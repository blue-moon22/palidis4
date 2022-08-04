import argparse, os
import unittest
from unittest.mock import patch, call, ANY
from pathlib import Path

from bin.get_IS_info import *

class TestGetISInfo(unittest.TestCase):
    TEST_COBS_TABLE = 'tests/data/input/test_insertion_sequences.fasta_1_results_table.txt'
    TEST_FFQ_JSON_LOC = 'tests/data/input'
    TEST_TAB_FILE = 'tests/data/input/test_insertion_sequence_annotations.tab'
    TEST_OUTPUT_PREFIX = 'tests/data/output/test'

    def test_get_cobs_info(self):
        actual = get_cobs_info(self.TEST_COBS_TABLE, self.TEST_FFQ_JSON_LOC)

        self.maxDiff = None
        self.assertEqual(actual,
        {'IS_cluster_115105_length_655': [('SAMN07658209', 'Mycobacterium tuberculosis')],
        'IS_cluster_2614_length_1455': [('SAMN07658618', 'Mycobacterium tuberculosis')]})

    def test_write_info(self):
        cobs_info = get_cobs_info(self.TEST_COBS_TABLE, self.TEST_FFQ_JSON_LOC)

        actual = write_info(self.TEST_TAB_FILE, cobs_info, self.TEST_OUTPUT_PREFIX)

        tab_name = self.TEST_OUTPUT_PREFIX + '_insertion_sequences_info.txt'
        tab = open(tab_name, "r")
        actual = "".join(tab.readlines())
        os.remove(tab_name)
        self.maxDiff = None
        self.assertEqual(actual, """IS_name\tsample_id\tcontig\titr1_start_position\titr1_end_position\titr2_start_position\titr2_end_position\titr_cluster\tCOBS_index_biosample_id\tCOBS_index_origin\nIS_cluster_115105_length_655\tSRS013170\tNODE_18_length_76504_cov_9.77495\t74408\t74436\t75032\t75062\t115105\tSAMN07658209\tMycobacterium tuberculosis\nIS_cluster_2614_length_1455\tSRS013170\tNODE_31_length_64375_cov_7.58579\t10034\t10063\t11459\t11488\t2614\tSAMN07658618\tMycobacterium tuberculosis\nIS_cluster_19090_length_1430\tSRS013170\tNODE_51_length_53875_cov_10.3394\t16273\t16301\t17674\t17702\t19090\t\t\n""")

    def test_arguments(self):
        actual = get_arguments().parse_args(
            ['--tab_file', 'tab_file', '--cobs_search_out', 'cobs_search_out',
            '--ffq_json', 'ffq_json', '--output_prefix', 'output_prefix'])
        self.assertEqual(actual, argparse.Namespace(tab_file='tab_file',
        ffq_json='ffq_json', cobs_table='cobs_search_out',
        output_prefix='output_prefix'))

    @patch('bin.get_IS_info.get_cobs_info')
    @patch('bin.get_IS_info.write_info')
    def test_main(self, mock_write_info, mock_get_cobs_info):
        args = get_arguments().parse_args(
            ['--tab_file', 'tab_file',
            '--cobs_search_out', 'cobs_search_out',
            '--ffq_json', 'ffq_json', '--output_prefix', 'output_prefix'])
        mock_get_cobs_info.return_value = {'final dictionary'}

        main(args)

        self.assertEqual(mock_get_cobs_info.call_args_list, [call('cobs_search_out', 'ffq_json')])
        self.assertEqual(mock_write_info.call_args_list, [call('tab_file', {'final dictionary'}, 'output_prefix')])
