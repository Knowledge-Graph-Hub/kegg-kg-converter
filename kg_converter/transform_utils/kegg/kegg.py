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
        source_name = "kegg"
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

        # Tables
        cpd_list = 'compounds.tsv'
        rn_list = 'reactions.tsv'
        path_list = 'pathways.tsv'
        ko_list = 'ko.tsv'
        cpd_path_link = 'compoundPathwayLink.tsv'
        cpd_rn_link = 'compoundReactionLink.tsv'
        rn_path_link = 'reactionPathwayLink.tsv'
        ko_path_link = 'koPathwayLink.tsv'
        ko_rn_link = 'koReactionLink.tsv'
        cpd2chebi = 'cpd2chebi.tsv'
        full_cpd= 'kegg-compounds.tsv'
        full_rn = 'kegg-reactions.tsv'
        full_path = 'kegg-pathways.tsv'

        

        # Pandas DF of 'list' files
        cpd_list_df = pd.read_csv(os.path.join(self.input_base_dir,cpd_list), low_memory=False, sep='\t')
        path_list_df = pd.read_csv(os.path.join(self.input_base_dir,path_list), low_memory=False, sep='\t')
        rn_list_df = pd.read_csv(os.path.join(self.input_base_dir, rn_list), low_memory=False, sep='\t')
        ko_list_df = pd.read_csv(os.path.join(self.input_base_dir, ko_list), low_memory=False, sep='\t')

        # transform data, something like:
        with open(os.path.join(self.input_base_dir,cpd_path_link), 'r') as cplf, \
                open(os.path.join(self.input_base_dir,rn_path_link), 'r') as rplf, \
                open(os.path.join(self.input_base_dir,cpd_rn_link), 'r') as crlf, \
                open(self.output_node_file, 'w') as node, \
                open(self.output_edge_file, 'w') as edge:

                # write headers (change default node/edge headers if necessary
                node.write("\t".join(self.node_header) + "\n")
                edge.write("\t".join(self.edge_header) + "\n")
                
                header_items = parse_header(cplf.readline(), sep='\t')
                
                seen_node: dict = defaultdict(int)
                seen_edge: dict = defaultdict(int)

                # Nodes
                cpd_node_type = "biolink:ChemicalSubstance"
                path_node_type = "biolink:Pathway"
                rn_node_type = "biolink:MolecularActivity"
                ko_node_type = "biolink:GeneFamily"

                # Edges
                cpd_to_path_label = "biolink:ChemicalToPathwayAssociation"
                cpd_to_rn_label = "#_ChemicalToReactionAssociaton"
                path_to_rn_label = "#_PathToReactionAssociation"

                # Compound & Pathway nodes
                for line in cplf:
                    # transform line into nodes and edges
                    # node.write(this_node1)
                    # node.write(this_node2)
                    # edge.write(this_edge)
                    items_dict = parse_line(line, header_items, sep='\t')
                    cpd_id = items_dict['cpdId']
                    cpd_name = cpd_list_df[cpd_list_df['cpdId'] == cpd_id]['cpd'].values[0]
                    
                    path_id = items_dict['pathwayId']
                    path_name = path_list_df[path_list_df['pathwayId'] == path_id]['pathway'].values[0]

                    if cpd_id not in seen_node:
                        write_node_edge_item(fh=node,
                                                header=self.node_header,
                                                data=[cpd_id,
                                                      cpd_name,
                                                      cpd_node_type])
                        seen_node[cpd_id] += 1


                    if path_id not in seen_node:
                        write_node_edge_item(fh=node,
                                                header=self.node_header,
                                                data=[path_id,
                                                        path_name,
                                                        path_node_type])
                        seen_node[path_id] += 1

                    #import pdb; pdb.set_trace()

        return None




        





        