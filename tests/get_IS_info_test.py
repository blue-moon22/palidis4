import argparse, os
import unittest
from unittest.mock import patch, call, ANY
from pathlib import Path

from bin.get_IS_info import *

class TestGetISInfo(unittest.TestCase):
    TEST_BLAST_OUT = 'tests/data/input/test_blast.out'
    TEST_ISFINDER_INFO = 'tests/data/input/IS_test.csv'
    TEST_COBS_TABLE = 'tests/data/input/test_insertion_sequences.fasta_1_results_table.txt'
    TEST_TAB_FILE = 'tests/data/input/test_insertion_sequence_annotations.tab'
    TEST_OUTPUT_PREFIX = 'tests/data/output/test'

    def test_create_annotations_dict(self):
        actual = create_annotations_dict(self.TEST_BLAST_OUT)

        self.maxDiff = None
        self.assertEqual(actual,
        {'IS_cluster_115105_length_655': [('IS_family', 'IS256')],
        'IS_cluster_2614_length_1455': [('IS', 'IS1249')]})

    def test_get_isfinder_info(self):
        actual = get_isfinder_info(self.TEST_ISFINDER_INFO)

        self.assertEqual(actual['IS1249'], 'Corynebacterium xerosis')

    def test_get_cobs_info(self):
        actual = get_cobs_info(self.TEST_COBS_TABLE)

        self.maxDiff = None
        self.assertEqual(actual,
        {'IS_cluster_2614_length_1455':
            [('SAMN09914449', 'Streptococcus mitis'),
            ('SAMN03956209', 'Streptococcus mitis'),
            ('SAMEA104055241', 'Streptococcus viridans'),
            ('SAMEA4555960', 'Gardnerella vaginalis'),
            ('SAMEA3319306', 'Streptococcus pneumoniae'),
            ('SAMN07658190', 'Mycobacterium tuberculosis')],
        'IS_cluster_115105_length_655':
            [('SAMEA709093', 'Haemophilus parainfluenzae'),
            ('SAMN09011143', 'Haemophilus haemolyticus')]})

    def test_write_info(self):
        is_finder_annot = create_annotations_dict(self.TEST_BLAST_OUT)
        is_finder_info = get_isfinder_info(self.TEST_ISFINDER_INFO)
        cobs_info = get_cobs_info(self.TEST_COBS_TABLE)

        actual = write_info(self.TEST_TAB_FILE, is_finder_annot, is_finder_info, cobs_info, self.TEST_OUTPUT_PREFIX)

        tab_name = self.TEST_OUTPUT_PREFIX + '_insertion_sequences_info.txt'
        tab = open(tab_name, "r")
        actual = "".join(tab.readlines())
        os.remove(tab_name)
        self.assertEqual(actual, """IS_name\tsample_id\tcontig\titr1_start_position	itr1_end_position\titr2_start_position\titr2_end_position\titr_cluster\tISfinder_name\tISfinder_origin\tpredicted_IS_family\tCOB_index_biosample_id\tCOB_index_origin\nIS_cluster_115105_length_655\tSRS013170\tNODE_18_length_76504_cov_9.77495\t74408\t74436\t75032\t75062\t115105\t\t\tIS256\tSAMEA709093;SAMN09011143\tHaemophilus parainfluenzae;Haemophilus haemolyticus\nIS_cluster_2614_length_1455\tSRS013170\tNODE_31_length_64375_cov_7.58579\t10034\t10063\t11459\t11488\t2614\tIS1249\tCorynebacterium xerosis\t\tSAMN09914449;SAMN03956209;SAMEA104055241;SAMEA4555960;SAMEA3319306;SAMN07658190\tStreptococcus mitis;Streptococcus mitis;Streptococcus viridans;Gardnerella vaginalis;Streptococcus pneumoniae;Mycobacterium tuberculosis\nIS_cluster_19090_length_1430\tSRS013170\tNODE_51_length_53875_cov_10.3394\t16273\t16301\t17674\t17702\t19090\t\t\t\t\t\n""")

    def test_arguments(self):
        actual = get_arguments().parse_args(
            ['--blast_out', 'blast_out', '--tab_file', 'tab_file',
            '--is_finder_info', 'is_finder_info', '--cobs_search_out', 'cobs_search_out',
            '--output_prefix', 'output_prefix'])
        self.assertEqual(actual, argparse.Namespace(blast_out='blast_out', tab_file='tab_file',
        is_info_csv='is_finder_info', cobs_table='cobs_search_out',
        output_prefix='output_prefix'))

    @patch('bin.get_IS_info.create_annotations_dict')
    @patch('bin.get_IS_info.get_isfinder_info')
    @patch('bin.get_IS_info.get_cobs_info')
    @patch('bin.get_IS_info.write_info')
    def test_main(self, mock_write_info, mock_get_cobs_info, mock_get_isfinder_info, mock_create_annotations_dict):
        args = get_arguments().parse_args(
            ['--blast_out', 'blast_out', '--tab_file', 'tab_file',
            '--is_finder_info', 'is_finder_info', '--cobs_search_out', 'cobs_search_out',
            '--output_prefix', 'output_prefix'])
        mock_create_annotations_dict.return_value = {'a dictionary'}
        mock_get_isfinder_info.return_value = {'another dictionary'}
        mock_get_cobs_info.return_value = {'final dictionary'}

        main(args)

        self.assertEqual(mock_create_annotations_dict.call_args_list, [call('blast_out')])
        self.assertEqual(mock_get_isfinder_info.call_args_list, [call('is_finder_info')])
        self.assertEqual(mock_get_cobs_info.call_args_list, [call('cobs_search_out')])
        self.assertEqual(mock_write_info.call_args_list, [call('tab_file', {'a dictionary'}, {'another dictionary'}, {'final dictionary'}, 'output_prefix')])
