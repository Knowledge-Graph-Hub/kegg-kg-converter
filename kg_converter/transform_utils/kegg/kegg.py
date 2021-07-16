#import csvimport 
import os
from typing import Dict, List, Optional
from collections import defaultdict

from numpy.lib.shape_base import column_stack

from kg_converter.transform_utils.transform import Transform
from kg_converter.utils.transform_utils import parse_header, parse_line, write_node_edge_item

import pandas as pd

"""
Ingest data from the KEGG database

Data of interest:
    1.  Compounds
    2.  Reactions
    3.  Pathways
    4.  KEGG Orthology

"""

class KEGGTransform(Transform):

    def __init__(self, input_dir: str = None, output_dir: str = None, nlp = False) -> None:
        source_name = 'kegg'
        super().__init__(source_name, input_dir, output_dir, nlp)  # set some variables

        self.node_header = ['id', 'name', 'category', 'exact_match', 'close_match', 'description' ]
        self.edge_header = ['subject', 'predicate', 'object', 'relation']
        self.nlp = nlp
    
    def run(self, data_file: Optional[str] = None):
        """Method is called and performs needed transformations to process the 
        KEGG data, additional information on this data can be found in the comment 
        at the top of this script.

        :param input_dir: SOurce of the downloaded data from the donwload step.
        :return: None
        
        """

        

        # Pandas DF of 'list' files
        #cpd_list_df = pd.read_csv(self.cpd_list, low_memory=False, sep='\t')
        #path_list_df = pd.read_csv(self.path_list, low_memory=False, sep='\t')
        #rn_list_df = pd.read_csv(self.rn_list, low_memory=False, sep='\t')
        #ko_list_df = pd.read_csv(self.ko_list, low_memory=False, sep='\t')
        #cpd_to_chebi_df = pd.read_csv(self.cpd2chebi, low_memory=False, sep='\t')

        # Pandas DF of 'kegg-*.tsv' files
        
        path_df = self.prune_columns(pd.read_csv(self.full_path, low_memory=False, sep='\t', usecols=['ENTRY', 'NAME']), 'path')
        rn_df = self.prune_columns(pd.read_csv(self.full_rn, low_memory=False, sep='\t', usecols=['ENTRY', 'DEFINITION', 'EQUATION']), 'rn')
        ko_df = self.prune_columns(pd.read_csv(self.full_ko, low_memory=False, sep='\t', usecols=['ENTRY', 'DEFINITION', 'DBLINKS']), 'ko')
        
        ## **********************************************************************
        # Establishing 1-to-1 relation between KO and XRefs (['DBLINKS'] column)
        ##***********************************************************************

        # Explode DBLINKS in ko_df to separate rows
        ko_df['DBLINKS'] = ko_df['DBLINKS'].apply(lambda row : str(row).split('|'))
        ko_df = ko_df.explode('DBLINKS')

        #ko_df['ID'] = ko_df['ID'].apply(lambda row : 'ko:'+str(row))
        ko_df['DBLINKS'] = ko_df['DBLINKS'].apply(lambda row : str(row).replace('RN: ', 'KEGG.REACTION:'))
        ko_df['DBLINKS'] = ko_df['DBLINKS'].apply(lambda row : str(row).strip().replace('COG: ', 'COG:'))
        ko_df['DBLINKS'] = ko_df['DBLINKS'].apply(lambda row : str(row).strip().replace('GO: ', 'go:'))
        ko_df['DBLINKS'] = ko_df['DBLINKS'].apply(lambda row : str(row).strip().replace('TC: ', 'tcdb:'))
        ko_df['DBLINKS'] = ko_df['DBLINKS'].apply(lambda row : str(row).strip().replace('CAZy: ', 'cazy:'))
        ko_df['DBLINKS'] = ko_df['DBLINKS'].apply(lambda row : str(row).strip().replace('UniProt: ', 'uniprot:'))
        ko_df['DBLINKS'] = ko_df['DBLINKS'].apply(lambda row: str(row).split(' '))
        # Add prefixes to all DBLINKS
        ko_df['DBLINKS'] = ko_df['DBLINKS'] \
                            .apply(lambda row: [str(row[0])]+[str(row[0])
                                .split(':')[0] + ':'+ x \
                                    for x in row \
                                        if not str(x).startswith(str(row[0]).split(':')[0]+ ':')])

        ko_df['DBLINKS'] = ['|'.join(map(str, l)) for l in ko_df['DBLINKS']]
        # Roll up to consolidated rows
        ko_df = ko_df.groupby(['ID', 'DESCRIPTION'], as_index=False).agg({'DBLINKS': lambda x: '|'.join(x)})
        ##########################################################################

        node_dict: dict = defaultdict(int)
        edge_dict: dict = defaultdict(int)

        df_dict = {
            'pathway': path_df,
            'rn': rn_df,
            'ko': ko_df
        }

        node_dict, edge_dict = self.post_data(self.path_cpd_link, node_dict, edge_dict, df_dict, 'w')
        node_dict, edge_dict = self.post_data(self.rn_cpd_link, node_dict, edge_dict, df_dict, 'a')
        node_dict, edge_dict = self.post_data(self.path_rn_link, node_dict, edge_dict, df_dict, 'a')
        node_dict, edge_dict = self.post_data(self.path_ko_link, node_dict, edge_dict, df_dict, 'a')
        node_dict, edge_dict = self.post_data(self.rn_ko_link, node_dict, edge_dict, df_dict, 'a')
                    

        return None

    def post_data(self, file, seen_node, seen_edge, desc_dict, mode):
        '''
        This function transforms the following KEGG data into nodes and edges:
            -   Pathway <-> Compound
            -   Reaction <-> Compound
            -   Pathways <-> Reaction
            -   Pathway <-> KEGG Orthology
            -   Reaction <-> KEGG Orthology

        :param file: The link file used as input.
        :param seen_node: Dictionary of all nodes recorded to avoid duplication.
        :param seen_edge: Dictionary of all edges recorded to avoid duplication.
        :param mode: Two options ['write' and 'append'] to avoid overwriting of nodes and edges tsv files.
        :return: seen_node and seen_edge such that the nodes and edges are unique throughout the process.
        '''

        with open(file, 'r') as f, \
                open(self.output_node_file, mode) as node, \
                open(self.output_edge_file, mode) as edge:

                # write headers (change default node/edge headers if necessary
                if mode == 'w':
                    node.write('\t'.join(self.node_header) + '\n')
                    edge.write('\t'.join(self.edge_header) + '\n')
                
                seen_node: dict = defaultdict(int)
                seen_edge: dict = defaultdict(int)

                # Nodes
                cpd_node_type = 'biolink:ChemicalSubstance'
                path_node_type = 'biolink:Pathway'
                rn_node_type = 'biolink:MolecularActivity'
                ko_node_type = 'biolink:GeneFamily'

                # Node Prefixes
                cpd_pref = 'KEGG.COMPOUND:'
                rn_pref = 'KEGG.REACTION:'
                path_pref = 'KEGG.PATHWAY:'
                ko_pref = 'KEGG.ORTHOLOGY:'

                # Edges
                path_to_cpd_label = 'biolink:has_participant'
                rn_to_cpd_label = 'biolink:has_participant'
                path_to_rn_label = 'biolink:has_participant'
                path_to_ko_label = 'biolink:has_participant'
                rn_to_ko_label = 'biolink:has_participant'
                predicate = ''
                predicate_curie = ''
                

                path_to_cpd_relation = 'RO:0000057'
                rn_to_cpd_relation = 'RO:0000057'
                path_to_rn_relation = 'RO:0000057'
                path_to_ko_relation = 'RO:0000057'
                rn_to_ko_relation = 'RO:0000057'

                #cpd_to_chebi_df = pd.DataFrame()
                node_id = ''
                node_pref = ''
                synonyms = ''
                xrefs = ''
                

                header_items = parse_header(f.readline(), sep='\t')

                if all(x in header_items for x in ['cpdId', 'pathwayId']):
                    predicate = path_to_cpd_label
                    predicate_curie = path_to_cpd_relation
                elif all(x in header_items for x in ['cpdId', 'rnId']):
                    predicate = rn_to_cpd_label
                    predicate_curie = rn_to_cpd_relation
                elif all(x in header_items for x in ['pathwayId', 'rnId']):
                    predicate = path_to_rn_label
                    predicate_curie = path_to_rn_relation
                elif all(x in header_items for x in ['koId', 'pathwayId']):
                    predicate = path_to_ko_label
                    predicate_curie = path_to_ko_relation
                elif all(x in header_items for x in ['koId', 'rnId']):
                    predicate = rn_to_ko_label
                    predicate_curie = rn_to_ko_relation
                else:
                    print('Unexpected column names provided.')
                
                for line in f:
                    # transform line into nodes and edges
                    # node.write(this_node1)
                    # node.write(this_node2)
                    # edge.write(this_edge)
                    items_dict = parse_line(line, header_items, sep='\t')
                    
                    edge_id = ''
                    subject = ''
                    object = ''
                    
                    
                    for key in items_dict.keys():
                        list_df = pd.DataFrame()
                        node_type = ''
                        desc_df = pd.DataFrame()
                        description = ''
                        #need_chebi = False
                        if key[:-2] == 'cpd':
                            list_df = pd.read_csv(self.cpd_list, sep='\t', low_memory=False)
                            node_type = cpd_node_type
                            node_pref = cpd_pref
                            #cpd_to_chebi_df = pd.read_csv(self.cpd2chebi, low_memory=False, sep='\t')
                            need_chebi = True
                        elif key[:-2] == 'rn':
                            list_df = pd.read_csv(self.rn_list, sep='\t', low_memory=False)
                            node_type = rn_node_type
                            node_pref = rn_pref
                        elif key[:-2] == 'pathway':
                            list_df = pd.read_csv(self.path_list, sep='\t', low_memory=False)
                            node_type = path_node_type
                            node_pref = path_pref
                        elif key[:-2] == 'ko':
                            list_df = pd.read_csv(self.ko_list, sep='\t', low_memory=False)
                            node_type = ko_node_type
                            node_pref = ko_pref

                        # Get CHEBI equivalent if possible
                        #if need_chebi and (any(cpd_to_chebi_df[key].str.contains(items_dict[key]))):
                            #node_id = cpd_to_chebi_df[cpd_to_chebi_df[key]==items_dict[key]]['chebiId'].values[0]
                        #else:
                        core_id = items_dict[key].split(':')[1] #This is the id without any prefix
                        desc_df = desc_dict.get(key[:-2]) # get the relevant DataFrame for description
                        node_id = node_pref+core_id

                        # Get description for the node
                        if desc_df is not None and core_id in desc_df['ID'].values:
                            description = desc_df.loc[desc_df['ID'] == core_id, 'DESCRIPTION'].iloc[0]


                        if edge_id == '':
                            edge_id = node_id
                            subject = node_id
                        else:
                            edge_id += '-'+node_id
                            object = node_id
                        
                        

                        if key =='pathwayId':
                            if 'rn' in items_dict[key]:
                                names = list_df[list_df[key] == items_dict[key].replace('rn', 'map')][key[:-2]].values[0]
                            elif 'ko' in items_dict[key]:
                                names = list_df[list_df[key] == items_dict[key].replace('ko', 'map')][key[:-2]].values[0]
                            else:
                                names = list_df[list_df[key] == items_dict[key]][key[:-2]].values[0]

                        else:
                            names = list_df[list_df[key] == items_dict[key]][key[:-2]].values[0]
                            if key == 'koId':
                                
                                xrefs = desc_dict[key[:-2]][desc_dict[key[:-2]]['ID'] == core_id]['DBLINKS'].values[0]
                                
                        
                        name = names.split(';')[0]
                        synonyms = ' | '.join(names.split(';')[1:]).strip()
                        
                        

                        # Nodes
                        if node_id not in seen_node:
                            write_node_edge_item(fh=node,
                                                header=self.node_header,
                                                data=[node_id,
                                                      name,
                                                      node_type,
                                                      synonyms,
                                                      xrefs,
                                                      description])
                            seen_node[node_id] += 1
                        

                    # Edges
                    if edge_id not in seen_edge:
                        write_node_edge_item(fh=edge,
                                        header=self.edge_header,
                                        data=[subject,
                                            predicate,
                                            object,
                                            predicate_curie])
                        seen_edge[edge_id] += 1

        return [seen_node, seen_edge]


    def prune_columns(self, df:pd.DataFrame, type:str)->pd.DataFrame:
        column_names = ['ID', 'DESCRIPTION']
        new_df = pd.DataFrame(columns=column_names)
        new_df['ID'] = df['ENTRY'].str.split(' ').str[0]
        if type == 'path':
            new_df['DESCRIPTION'] = df['NAME'].str.split('DESCRIPTION').str[1]
        elif type == 'rn':
            new_df['DESCRIPTION'] = df['DEFINITION'] + ' | EQUATION: ' + df['EQUATION']
        elif type == 'ko':
            new_df['DESCRIPTION'] = df['DEFINITION']
            new_df['DBLINKS'] = df['DBLINKS']
        else:
            print('Unknown type of data')

        return new_df.dropna()



        





        