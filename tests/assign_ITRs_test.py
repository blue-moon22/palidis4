import argparse, os
import unittest
from unittest.mock import patch, call, ANY

from bin.assign_ITRs import create_cluster_dictionary, find_itr_clusters, write_fasta, get_arguments, main

class TestAssignITRs(unittest.TestCase):
    TEST_CLIPPED_READS = 'tests/data/input/test_irs.fasta'
    TEST_CLSTR1 = 'tests/data/input/test1.fasta.clstr'
    TEST_CLSTR2 = 'tests/data/input/test2.fasta.clstr'
    TEST_INFO_TAB1 = 'tests/data/input/test1_contigs_reads_ir_position_info.tab'
    TEST_INFO_TAB2 = 'tests/data/input/test2_contigs_reads_ir_position_info.tab'
    TEST_OUTPUT_PREFIX = 'tests/data/output/test'

    def test_create_cluster_dictionary(self):
        actual = create_cluster_dictionary(self.TEST_CLSTR1)
        self.assertEqual(actual[1][95]['ERR589346'], '0')
        self.assertEqual(actual[0][95]['ERR589347'], '1')
        self.assertEqual(actual[0][157]['ERR589346'], '0')

    def test_create_multi_cluster_dictionary(self):
        actual = create_cluster_dictionary(self.TEST_CLSTR2)
        self.assertEqual(actual[1][95]['ERR589346'], '0')
        self.assertEqual(actual[0][95]['ERR589347'], '1')
        self.assertEqual(actual[0][157]['ERR589346'], '0')
        self.assertEqual(actual[1][99]['ERR589347'], '1')

    def test_find_itr_clusters(self):
        cl_dict = create_cluster_dictionary(self.TEST_CLSTR1)
        find_itr_clusters(cl_dict, self.TEST_INFO_TAB1, self.TEST_OUTPUT_PREFIX)
        tab_name = self.TEST_OUTPUT_PREFIX + '_contigs_reads_itr_position_info.tab'
        tab = open(tab_name, "r")
        actual = "".join(tab.readlines())
        #os.remove(tab_name)
        self.assertEqual(actual, """sample_id\tcontig\tread1\tread2\tposition1\tposition2\titr_cluster\nERR589346\tNODE_3_length_222578_cov_8.12661\tSeq158_nstart_ERR589346_nend_ERR589346.133_FCC4C01ACXX:6:1101:12524:2220#AAGTCTCT/2_f2\tSeq96_nstart_ERR589346_nend_ERR589346.126_FCC4C01ACXX:6:1101:12254:2132#AAGTCTCT/2_f1\t100\t600\t0\n""")


    def test_find_multi_itr_clusters(self):
        cl_dict = create_cluster_dictionary(self.TEST_CLSTR2)
        find_itr_clusters(cl_dict, self.TEST_INFO_TAB2, self.TEST_OUTPUT_PREFIX + '_multi')
        tab_name = self.TEST_OUTPUT_PREFIX + '_multi_contigs_reads_itr_position_info.tab'
        tab = open(tab_name, "r")
        actual = "".join(tab.readlines())
        #os.remove(tab_name)
        self.assertEqual(actual, """sample_id\tcontig\tread1\tread2\tposition1\tposition2\titr_cluster\nERR589346\tNODE_3_length_222578_cov_8.12661\tSeq158_nstart_ERR589346_nend_ERR589346.133_FCC4C01ACXX:6:1101:12524:2220#AAGTCTCT/2_f2\tSeq96_nstart_ERR589346_nend_ERR589346.126_FCC4C01ACXX:6:1101:12254:2132#AAGTCTCT/2_f1\t100\t600\t0\nERR589347\tNODE_1_length_300000_cov_1.10001\tSeq96_nstart_ERR589347_nend_ERR589347.195_FCC4C01ACXX:6:1101:18660:2181#AAGTCTCT/1_f2\tSeq100_nstart_ERR589347_nend_ERR589347.195_FCC4C01ACXX:6:1101:18660:2181#AAGTCTCT/2_f1\t500\t1100\t1\n""")


    def test_arguments(self):
        actual = get_arguments().parse_args(
            ['--clipped_reads', 'clipped_reads', '--cdhit_cluster_file', 'clstr_file', '--info_tab_file', 'tab_file',
            '--output_prefix', 'output_prefix'])
        self.assertEqual(actual,
                         argparse.Namespace(clipped_fasta='clipped_reads', cluster_file='clstr_file', tab_file='tab_file', output_prefix='output_prefix'))

    def test_arguments_short_options(self):
        actual = get_arguments().parse_args(
            ['-r', 'clipped_reads', '-c', 'clstr_file', '-t', 'tab_file',
            '-o', 'output_prefix'])
        self.assertEqual(actual,
                         argparse.Namespace(clipped_fasta='clipped_reads', cluster_file='clstr_file', tab_file='tab_file', output_prefix='output_prefix'))

    def test_write_fasta(self):
        alloc_1 = [0]*100
        alloc_1[98] = 1
        alloc_2 = [0]*100
        alloc_2[22] = 1
        write_fasta(self.TEST_CLIPPED_READS, (alloc_1, alloc_1), self.TEST_OUTPUT_PREFIX)
        fasta_name = self.TEST_OUTPUT_PREFIX + '_ITRs.fasta'
        tab = open(fasta_name, "r")
        actual = "".join(tab.readlines())
        os.remove(fasta_name)
        self.assertEqual(actual, """>Seq99_nstart_ERR589535_nend_ERR589535.33_f1_LCoord_29_RCoord_56\nGAGATCGGCCCGGCCTTTGAGCACCACG\n""")
        write_fasta(self.TEST_CLIPPED_READS, (alloc_1, alloc_2), self.TEST_OUTPUT_PREFIX)

    @patch('bin.assign_ITRs.create_cluster_dictionary')
    @patch('bin.assign_ITRs.find_itr_clusters')
    @patch('bin.assign_ITRs.write_fasta')
    def test_main(self, mock_write_fasta, mock_find_itr_clusters, mock_create_cluster_dictionary):
        args = get_arguments().parse_args(
            ['--clipped_reads', 'clipped_fasta', '--cdhit_cluster_file', 'clstr_file', '--info_tab_file', 'tab_file',
            '--output_prefix', 'output_prefix'])
        mock_create_cluster_dictionary.return_value = [{}]
        mock_find_itr_clusters.return_value = {}

        main(args)

        self.assertEqual(mock_create_cluster_dictionary.call_args_list, [call('clstr_file')])
        self.assertEqual(mock_find_itr_clusters.call_args_list, [call([{}], 'tab_file', 'output_prefix')])
        self.assertEqual(mock_write_fasta.call_args_list, [call('clipped_fasta', {}, 'output_prefix')])
