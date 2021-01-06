import csvimport os
from typing import Dict, List, Optional
from collections import defaultdict

from kg_converter.transform_utils.transform import Transform
from kg_converter.utils.transform_utils import parse_header, parse_line, write_node_edge_item

"""
Ingest data from the KEGG database

Data of interest:
    1.  Compounds
    2.  Reactions
    3.  Pathways

"""

class KEGGTransform(Transform):

    def __init__(self, input_dir: str = None, output_dir: str = None, nlp = False) -> None:
        source_name = "kegg"
        super().__init__(source_name, input_dir, output_dir, nlp)  # set some variables

        self.node_header = ['id', 'name', 'category', 'curie']
        self.edge_header = ['subject', 'edge_label', 'object', 'relation']
        self.nlp = nlp

    def run(self, data_file: Optional[str] = None):
        """Method is called and performs needed transformations to process the 
        KEGG data, additional information on this data can be found in the comment 
        at the top of this script"""

        