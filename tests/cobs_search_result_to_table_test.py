import argparse, os
import unittest
from unittest.mock import patch, call, ANY

from bin.cobs_search_result_to_table import *

class TestCobsSearchResultToTable(unittest.TestCase):
    TEST_FAI = "tests/data/input/test_insertion_sequences.fasta.fai"
    TEST_COBS_RESULTS = "tests/data/input/test_insertion_sequences.fasta_1_results.txt"
    TEST_OUTPUT = "tests/data/output/test_insertion_sequences.fasta_1_results_table.txt"
    def test_get_query_size(self):
        actual = get_query_size(self.TEST_FAI, 1)

        self.maxDiff = None
        self.assertEqual(actual['IS_cluster_107290_length_2079'], {'hits': {}, 'max_kmers': 2079})

    def test_get_results(self):
        queries = get_query_size(self.TEST_FAI, 1)

        actual = get_results(queries, self.TEST_COBS_RESULTS)

        self.assertEqual(actual['IS_cluster_107290_length_2079'], {'max_kmers': 2079, 'hits': {}})
        self.assertEqual(actual['IS_cluster_222669_length_1855'], {'max_kmers': 1855, 'hits': {'SAMN07658209': 1.0}})

    def test_write_table(self):
        queries = get_query_size(self.TEST_FAI, 1)
        results = get_results(queries, self.TEST_COBS_RESULTS)

        actual = write_table(results, self.TEST_OUTPUT)

        tab = open(self.TEST_OUTPUT, "r")
        actual = "".join(tab.readlines())
        os.remove(self.TEST_OUTPUT)
        self.assertEqual(actual, """query\tsample_id\tkmer_similarity\nIS_cluster_222669_length_1855\tSAMN07658209\t1.0\nIS_cluster_50504_length_2145\tSAMN07658618\t1.0\n""")

    def test_arguments(self):
        actual = get_arguments().parse_args(
            ['--cobs_outfile', 'cobs_outfile', '--fai_file', 'fai_file', '--outname', 'outname', '--kmer_length', '1'])
        self.assertEqual(actual, argparse.Namespace(cobs_outfile='cobs_outfile', fai_file='fai_file', outname='outname', kmer_length=1))

    @patch('bin.cobs_search_result_to_table.get_query_size')
    @patch('bin.cobs_search_result_to_table.get_results')
    @patch('bin.cobs_search_result_to_table.write_table')
    def test_main(self, mock_write_table, mock_get_results, mock_get_query_size):
        args = get_arguments().parse_args(
            ['--cobs_outfile', 'cobs_outfile', '--fai_file', 'fai_file', '--outname', 'outname', '--kmer_length', '1'])
        mock_get_query_size.return_value = {'a dictionary'}
        mock_get_results.return_value = {'another dictionary'}

        main(args)

        self.assertEqual(mock_get_query_size.call_args_list, [call('fai_file', 1)])
        self.assertEqual(mock_get_results.call_args_list, [call({'a dictionary'}, 'cobs_outfile')])
        self.assertEqual(mock_write_table.call_args_list, [call({'another dictionary'}, 'outname')])
