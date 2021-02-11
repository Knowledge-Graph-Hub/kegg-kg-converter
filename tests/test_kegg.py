import os
import unittest
import pandas as pd
from kg_converter.transform_utils.kegg import KEGGTransform
from kg_converter.transform_utils.kegg.kegg import parse_line
from kg_converter.utils.transform_utils import parse_header
from parameterized import parameterized


class TestKEGG(unittest.TestCase):

    def setUp(self) -> None:
        self.resources = 'tests/resources/'
        self.input_dir = os.path.join(self.resources,'kegg/input/')
        self.output_dir = os.path.join(self.resources,'kegg/output/')
        self.test_file = 'kegg_test.tsv'
        self.kegg_fh = open(os.path.join(self.input_dir, self.test_file),'rt')
        self.kegg_output_dir = os.path.join(self.output_dir, "kegg/")
        self.kegg = KEGGTransform(input_dir=self.input_dir,
                                                 output_dir=self.output_dir)

    @parameterized.expand([
     ('subjectId', 'path:map00010'),
     ('objectId', 'cpd:C00022'),
     ])
    def test_parse_kegg_line(self, key, value):
        header = parse_header(self.kegg_fh.readline())
        line = self.kegg_fh.readline()
        parsed = parse_line(line, header)
        self.assertTrue(key in parsed)
        self.assertEqual(value, parsed[key])

    def test_run(self):
        self.assertTrue(isinstance(self.kegg.run, object))
        self.kegg.run(data_file=self.test_file)
        self.assertTrue(os.path.isdir(self.kegg_output_dir))

    def test_nodes_file(self):
        self.kegg.run(data_file=self.test_file)
        node_file = os.path.join(self.kegg_output_dir, "nodes.tsv")
        self.assertTrue(os.path.isfile(node_file))
        node_df = pd.read_csv(node_file, sep="\t", header=0)
        self.assertEqual((23, 5), node_df.shape)
        self.assertEqual(['id', 'name', 'category', 'exact_match', 'description' ],
                         list(node_df.columns))

    def test_nodes_are_not_repeated(self):
        self.kegg.run(data_file=self.test_file)
        node_file = os.path.join(self.kegg_output_dir, "nodes.tsv")
        node_df = pd.read_csv(node_file, sep="\t", header=0)
        nodes = list(node_df.id)
        unique_nodes = list(set(nodes))
        self.assertCountEqual(nodes, unique_nodes)

    def test_edges_file(self):
        self.kegg.run(data_file=self.test_file)
        edge_file = os.path.join(self.dc_output_dir, "edges.tsv")
        self.assertTrue(os.path.isfile(edge_file))
        edge_df = pd.read_csv(edge_file, sep="\t", header=0)
        self.assertEqual((21, 7), edge_df.shape)
        self.assertEqual(
            ['subject', 'predicate', 'object', 'relation'],
             list(edge_df.columns)
        )
