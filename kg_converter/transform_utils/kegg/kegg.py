#import csvimport 
import os
from typing import Dict, List, Optional
from collections import defaultdict

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

        self.node_header = ['id', 'name', 'category']
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
        cpd_list_df = pd.read_csv(self.cpd_list, low_memory=False, sep='\t')
        path_list_df = pd.read_csv(self.path_list, low_memory=False, sep='\t')
        rn_list_df = pd.read_csv(self.rn_list, low_memory=False, sep='\t')
        ko_list_df = pd.read_csv(self.ko_list, low_memory=False, sep='\t')

        node_dict: dict = defaultdict(int)
        edge_dict: dict = defaultdict(int)

        node_dict, edge_dict = self.post_data(self.cpd_path_link, node_dict, edge_dict)
        #node_dict, edge_dict = self.post_data(self.cpd_rn_link, node_dict, edge_dict)
                    

        return None

    def post_data(self, file, seen_node, seen_edge):

        with open(file, 'r') as f, \
                open(self.output_node_file, 'w') as node, \
                open(self.output_edge_file, 'w') as edge:

                # write headers (change default node/edge headers if necessary
                node.write('\t'.join(self.node_header) + '\n')
                edge.write('\t'.join(self.edge_header) + '\n')
                
                seen_node: dict = defaultdict(int)
                seen_edge: dict = defaultdict(int)

                # Nodes
                cpd_node_type = 'biolink:ChemicalSubstance'
                path_node_type = 'biolink:Pathway'
                rn_node_type = 'biolink:MolecularActivity'
                ko_node_type = 'biolink:GeneFamily'

                # Edges
                cpd_to_path_label = 'biolink:ChemicalToPathwayAssociation'
                cpd_to_rn_label = '#_ChemicalToReactionAssociaton'
                path_to_rn_label = '#_PathToReactionAssociation'
                ko_to_path_label = ''
                ko_to_rn_label = ''
                predicate = ''
                predicate_curie = ''
                edge_id = ''
                subject = ''
                object = ''

                cpd_to_path_relation = 'NEED_CURIE'
                cpd_to_rn_relation = 'NEED_CURIE'
                path_to_rn_relation = 'NEED_CURIE'
                ko_to_path_relation = 'NEED_CURIE'
                ko_to_rn_relation = 'NEED_CURIE'

                header_items = parse_header(f.readline(), sep='\t')

                if all(x in header_items for x in ['cpdId', 'pathwayId']):
                    predicate = cpd_to_path_label
                    predicate_curie = cpd_to_path_relation
                elif all(x in header_items for x in ['cpdId', 'rnId']):
                    predicate = cpd_to_rn_label
                    predicate_curie = cpd_to_rn_relation
                elif all(x in header_items for x in ['pathwayId', 'rnId']):
                    predicate = path_to_rn_label
                    predicate_curie = path_to_rn_relation
                
                for line in f:
                    # transform line into nodes and edges
                    # node.write(this_node1)
                    # node.write(this_node2)
                    # edge.write(this_edge)
                    items_dict = parse_line(line, header_items, sep='\t')
                    
                    
                    for key in items_dict.keys():
                        list_df = pd.DataFrame()
                        node_type = ''
                        if key[:-2] == 'cpd':
                            list_df = pd.read_csv(self.cpd_list, sep='\t', low_memory=False)
                            node_type = cpd_node_type
                        elif key[:-2] == 'rn':
                            list_df = pd.read_csv(self.rn_list, sep='\t', low_memory=False)
                            node_type = rn_node_type
                        elif key[:-2] == 'pathway':
                            list_df = pd.read_csv(self.path_list, sep='\t', low_memory=False)
                            node_type = path_node_type
                        elif key[:-2] == 'ko':
                            list_df = pd.read_csv(self.ko_list, sep='\t', low_memory=False)
                            node_type = ko_node_type

                        node_id = items_dict[key]

                        if edge_id == '':
                            edge_id = node_id
                            subject = node_id
                        else:
                            edge_id += node_id
                            object = node_id

                        name = list_df[list_df[key] == node_id][key[:-2]].values[0]
                        # Nodes
                        if node_id not in seen_node:
                            write_node_edge_item(fh=node,
                                                header=self.node_header,
                                                data=[node_id,
                                                      name,
                                                      node_type])
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



        





        