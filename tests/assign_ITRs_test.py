import argparse, os
import unittest
from unittest.mock import patch, call, ANY

from bin.assign_ITRs import create_cluster_dictionary, find_itr_clusters, get_arguments, main

class TestAssignITRs(unittest.TestCase):
    TEST_FASTA = 'tests/data/input/test_IR_1.fasta'
    TEST_CLSTR1 = 'tests/data/input/test1.fasta.clstr'
    TEST_CLSTR2 = 'tests/data/input/test2.fasta.clstr'
    TEST_INFO_TAB1 = 'tests/data/input/test1_contigs_reads_ir_position_info.tab'
    TEST_INFO_TAB2 = 'tests/data/input/test2_contigs_reads_ir_position_info.tab'
    TEST_OUTPUT_PREFIX = 'tests/data/output/test'

    def test_create_cluster_dictionary(self):
        actual = create_cluster_dictionary(self.TEST_FASTA , self.TEST_CLSTR1)
        self.assertEqual(actual[1][95]['ERR589346'], '0')
        self.assertEqual(actual[0][95]['ERR589347'], '1')
        self.assertEqual(actual[0][157]['ERR589346'], '0')

    def test_create_multi_cluster_dictionary(self):
        actual = create_cluster_dictionary(self.TEST_FASTA , self.TEST_CLSTR2)
        self.assertEqual(actual[1][95]['ERR589346'], '0')
        self.assertEqual(actual[0][95]['ERR589347'], '1')
        self.assertEqual(actual[0][157]['ERR589346'], '0')
        self.assertEqual(actual[1][99]['ERR589347'], '1')

    def test_find_itr_clusters(self):
        cl_dict = create_cluster_dictionary(self.TEST_FASTA , self.TEST_CLSTR1)
        find_itr_clusters(cl_dict, self.TEST_INFO_TAB1, self.TEST_OUTPUT_PREFIX)
        tab_name = self.TEST_OUTPUT_PREFIX + '_contigs_reads_itr_position_info.tab'
        tab = open(tab_name, "r")
        actual = "".join(tab.readlines())
        os.remove(tab_name)
        self.assertEqual(actual, """sample_id\tcontig\tread1\tread2\tposition1\tposition2\titr_cluster\nERR589346\tNODE_3_length_222578_cov_8.12661\tSeq158_nstart_ERR589346_nend_ERR589346.133_FCC4C01ACXX:6:1101:12524:2220#AAGTCTCT/2_f2\tSeq96_nstart_ERR589346_nend_ERR589346.126_FCC4C01ACXX:6:1101:12254:2132#AAGTCTCT/2_f1\t100\t600\t0\n""")


    def test_find_multi_itr_clusters(self):
        cl_dict = create_cluster_dictionary(self.TEST_FASTA , self.TEST_CLSTR2)
        find_itr_clusters(cl_dict, self.TEST_INFO_TAB2, self.TEST_OUTPUT_PREFIX + '_multi')
        tab_name = self.TEST_OUTPUT_PREFIX + '_multi_contigs_reads_itr_position_info.tab'
        tab = open(tab_name, "r")
        actual = "".join(tab.readlines())
        os.remove(tab_name)
        self.assertEqual(actual, """sample_id\tcontig\tread1\tread2\tposition1\tposition2\titr_cluster\nERR589346\tNODE_3_length_222578_cov_8.12661\tSeq158_nstart_ERR589346_nend_ERR589346.133_FCC4C01ACXX:6:1101:12524:2220#AAGTCTCT/2_f2\tSeq96_nstart_ERR589346_nend_ERR589346.126_FCC4C01ACXX:6:1101:12254:2132#AAGTCTCT/2_f1\t100\t600\t0\nERR589347\tNODE_1_length_300000_cov_1.10001\tSeq96_nstart_ERR589347_nend_ERR589347.195_FCC4C01ACXX:6:1101:18660:2181#AAGTCTCT/1_f2\tSeq100_nstart_ERR589347_nend_ERR589347.195_FCC4C01ACXX:6:1101:18660:2181#AAGTCTCT/2_f1\t500\t1100\t1\n""")


    def test_arguments(self):
        actual = get_arguments().parse_args(
            ['--combined_itr_fasta', 'itr_fasta', '--cdhit_cluster_file', 'clstr_file', '--combined_info_tab_file', 'tab_file',
            '--output_prefix', 'output_prefix'])
        self.assertEqual(actual,
                         argparse.Namespace(itr_fasta='itr_fasta', cluster_file='clstr_file', tab_file='tab_file', output_prefix='output_prefix'))

    def test_arguments_short_options(self):
        actual = get_arguments().parse_args(
            ['-i', 'itr_fasta', '-c', 'clstr_file', '-t', 'tab_file',
            '-o', 'output_prefix'])
        self.assertEqual(actual,
                         argparse.Namespace(itr_fasta='itr_fasta', cluster_file='clstr_file', tab_file='tab_file', output_prefix='output_prefix'))

    @patch('bin.assign_ITRs.create_cluster_dictionary')
    @patch('bin.assign_ITRs.find_itr_clusters')
    def test_main(self, mock_find_itr_clusters, mock_create_cluster_dictionary):
        args = get_arguments().parse_args(
            ['--combined_itr_fasta', 'itr_fasta', '--cdhit_cluster_file', 'clstr_file', '--combined_info_tab_file', 'tab_file',
            '--output_prefix', 'output_prefix'])
        mock_create_cluster_dictionary.return_value = [{}]
        mock_find_itr_clusters.return_value = {}

        main(args)

        self.assertEqual(mock_create_cluster_dictionary.call_args_list, [call('itr_fasta', 'clstr_file')])
        self.assertEqual(mock_find_itr_clusters.call_args_list, [call([{}], 'tab_file', 'output_prefix')])
