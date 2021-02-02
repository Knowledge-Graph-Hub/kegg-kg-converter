#!/usr/bin/env python
# -*- coding: utf-8 -*-


import logging
import os
from urllib.request import Request, urlopen

import yaml
from os import path
from tqdm.auto import tqdm  # type: ignore
import pandas as pd
import urllib3
import io
import time

def parse_response(url: str)-> pd.DataFrame:
    '''
    Implement the KEGG API to get 'LIST', 'LINK' and 'CONV' of relevant information and returns a pandas dataframe.

    :param url: URL of the REST API
    :return: Pandas DataFrame.
    '''
    col_1 = ''
    col_2 = ''

    if url.split('/')[3] == 'list':
        col_1 = url.split('/')[4]+'Id'
        col_2 = url.split('/')[4]
    elif url.split('/')[3] == 'link':
        col_1 = url.split('/')[5]+'Id'
        col_2 = url.split('/')[4]+'Id'
    elif url.split('/')[3] == 'conv':
        col_1 = url.split('/')[5]+'Id'
        col_2 = url.split('/')[4]+'Id'
    
    cols = [col_1, col_2]
    df = pd.DataFrame(columns=cols)
    http = urllib3.PoolManager()
    pathway_response = http.request('GET', url, preload_content=False)
    pathway_response.auto_close=False

    for line in io.TextIOWrapper(pathway_response):
        df = df.append(pd.Series(line.strip('\n').split('\t'),index = df.columns), ignore_index=True)
    return df


def has_digit(string: str) -> bool:
    '''
    Simple function to check if a string is alphanumeric or not.

    :param string: String
    :return: bool (True/False)
    
    '''
    for s in string:
        if s.isdigit():
            return True
    return False

def parse_response_get(url: str, output_dir: str, fn: str)-> pd.DataFrame:
    '''
    Implement the KEGG API to 'GET' relevant information and returns a pandas dataframe.

    :param url: URL of the REST API
    :param output_dir: The target location for the data file
    :param fn: Filename of the target.

    :return: Pandas DataFrame

    '''
    list_of_dict = []
    non_column_chars = ['-', ' ', ';']
    file_path = os.path.join(output_dir, fn+'.tsv')
    df = pd.read_csv(file_path, sep='\t')
    df_id_col = df[df.columns[0]]

    
    for id in df_id_col.iteritems():
        dictionary = {}
        last_key = ''
        new_url = url+id[1]
        http = urllib3.PoolManager()
        pathway_response = http.request('GET', new_url, preload_content=False)
        pathway_response.auto_close=False
        
        for line in io.TextIOWrapper(pathway_response):
            line_elements = line.split('  ')
            list_of_elements = [x.strip() for x in line_elements if x]


            if list_of_elements[0].isupper() \
                and not has_digit(list_of_elements[0]) \
                and not any(map(list_of_elements[0].__contains__, non_column_chars)) \
                and len(list_of_elements) > 1 \
                and len(list_of_elements[0]) > 3:

                last_key = list_of_elements[0]
                        
                if last_key == 'ENZYME':
                    dictionary[last_key] = ' | '.join(list_of_elements[1:])
                elif last_key in dictionary.keys():
                    dictionary[last_key] += (' | '+'-'.join(list_of_elements[1:]))
                else:
                    dictionary[last_key] = ' '.join(list_of_elements[1:])

            else:
                if last_key == 'COMMENT':
                        dictionary[last_key] += (' '+' '.join(list_of_elements))
                else:
                    dictionary[last_key] += (' | '+'-'.join(list_of_elements))

            dictionary[last_key] = dictionary[last_key].replace(' | ///', '')

        list_of_dict.append(dictionary)

    return pd.DataFrame(list_of_dict)



def download_from_yaml(yaml_file: str, output_dir: str,
                       ignore_cache: bool = False) -> None:
    """Given an download info from an download.yaml file, download all files

    :param yaml_file: A string pointing to the download.yaml file, to be parsed for things to download.
    :param output_dir: A string pointing to where to write out downloaded files.
    :param ignore_cache: Ignore cache and download files even if they exist [false]

    :return: None.
    """

    os.makedirs(output_dir, exist_ok=True)
    filename_list = []
    nap_time = 10
    with open(yaml_file) as f:
        data = yaml.load(f, Loader=yaml.FullLoader)
        for item in tqdm(data, desc="Downloading files"):
            if 'url' not in item:
                logging.warning("Couldn't find url for source in {}".format(item))
                continue
            outfile = os.path.join(
                output_dir,
                item['local_name']
                if 'local_name' in item
                else item['url'].split("/")[-1]
            )
            logging.info("Retrieving %s from %s" % (outfile, item['url']))

            if path.exists(outfile):
                if ignore_cache:
                    logging.info("Deleting cached version of {}".format(outfile))
                    os.remove(outfile)
                else:
                    logging.info("Using cached version of {}".format(outfile))
                    continue
            
            url_breakdown = item['url'].split('/')
            data_name = '_'.join(url_breakdown[-2:])
            
            if url_breakdown[3] in ['list', 'link', 'conv']:
                # CONV
                if url_breakdown[3] == 'conv':
                    conv_list = [
                        ['cpd', 'chebi']
                    ]
                    
                    for c in conv_list:
                        new_url = item['url']+'/'.join(c)
                        fn = item['local_name'].replace('placeholder',('2').join(c))
                        if not path.exists(os.path.join(output_dir, fn)):
                            parse_response(new_url).to_csv(os.path.join(output_dir, fn), sep='\t', index=False)
                            # Uncomment below if len(conv_list) > 1
                            # print('Waiting a bit before next API call...')
                            # time.sleep(nap_time)
                # LIST or LINK
                else:
                    parse_response(item["url"]).to_csv(os.path.join(output_dir, item['local_name']), sep='\t', index=False)
                

            # GET        
            elif url_breakdown[3] in ['get']:
                list_of_dbs = ['pathways', 'reactions', 'compounds', 'ko']

                for element in list_of_dbs:
                    fn = item['local_name'].replace('placeholder', element)
                    print('Looking for '+fn+' ...')
                    if not path.exists(os.path.join(output_dir,fn)):
                        print('Not found. Getting: '+fn)
                        parse_response_get(item['url'], output_dir, element).to_csv(os.path.join(output_dir, fn), sep='\t', index=False)
                        print('Waiting a bit before next API call...')
                        time.sleep(nap_time)
                    else:
                        print('Found '+fn+'!')
            else:
                print('Non-KEGG URL')
                # OR regular wget download here.
                
                req = Request(item['url'], headers={'User-Agent': 'Mozilla/5.0'})
                with urlopen(req) as response, open(outfile, 'wb') as out_file:  # type: ignore
                        data = response.read()  # a `bytes` object
                        out_file.write(data)
                

    return None


