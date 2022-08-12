import argparse, os
import unittest
from unittest.mock import patch, call, ANY
from pathlib import Path

from bin.assign_ITRs import *

class TestAssignITRs(unittest.TestCase):
    TEST_CLIPPED_READS = 'tests/data/input/test_irs.fasta'
    TEST_CLSTR1 = 'tests/data/input/test.fasta.clstr'
    TEST_CLSTR2 = 'tests/data/input/test2.fasta.clstr'
    TEST_INFO_TAB1 = 'tests/data/input/test_contigs_reads_ir_position_info.tab'
    TEST_INFO_TAB2 = 'tests/data/input/test2_contigs_reads_ir_position_info.tab'
    TEST_OUTPUT_PREFIX = 'tests/data/output/test'
    TEST_CONTIG_FASTA = 'tests/data/input/test_assemblies_filtered.fasta'

    def test_assembly_bins(self):
        actual = create_assembly_bins(self.TEST_CONTIG_FASTA)

        self.assertEqual(len(actual['NODE_1_length_49561_cov_5.4235']), 49561)

    def test_create_cluster_dictionary(self):
        actual = create_cluster_dictionary(self.TEST_CLSTR1)
        self.assertEqual(actual[1][95], '')
        self.assertEqual(actual[0][95], '')
        self.assertEqual(actual[0][157], '')

    def test_bin_positions(self):
        assembly_bins_dict = create_assembly_bins(self.TEST_CONTIG_FASTA)
        cl_dict = create_cluster_dictionary(self.TEST_CLSTR1)
        actual = bin_positions(cl_dict, self.TEST_INFO_TAB1, assembly_bins_dict, self.TEST_OUTPUT_PREFIX)
        self.assertEqual(actual['NODE_823_length_1805_cov_1014.02']['0'][0][0], '0')

    def test_count_bins(self):
        actual = count_bins('0111110011110001')

        self.assertEqual(actual, [(2, 7), (9, 13), (16, 17)])

    def test_count_bins_with_Ns(self):
        actual = count_bins('100000000000011100')

        self.assertEqual(actual, [(1, 2), (14, 17)])

    def test_get_itr_sequences(self):
        assemblies_dict = create_assembly_bins(self.TEST_CONTIG_FASTA)

        actual = get_itr_sequences(assemblies_dict["NODE_823_length_1805_cov_1014.02"], (246, 271, 1242, 1266))

        self.assertEqual(actual, ['AATATTGTTTTACATCCTGCATCTTA', 'TAAGATGCAGGATGTAAAACAATAT'])

    @patch('bin.assign_ITRs.are_reverse_cmp')
    def test_write_itr_annotations(self, mock_are_reverse_cmp):
        mock_are_reverse_cmp.return_value = True
        assembly_bins_dict = create_assembly_bins(self.TEST_CONTIG_FASTA)
        cl_dict = create_cluster_dictionary(self.TEST_CLSTR1)
        clusters_positions = bin_positions(cl_dict, self.TEST_INFO_TAB1, assembly_bins_dict, self.TEST_OUTPUT_PREFIX)

        write_itr_annotations(clusters_positions, assembly_bins_dict, 500, 3000, 25, 50, self.TEST_OUTPUT_PREFIX, 8)

        tab_name = self.TEST_OUTPUT_PREFIX + '_insertion_sequence_annotations.tab'
        tab = open(tab_name, "r")
        actual = "".join(tab.readlines())
        os.remove(tab_name)
        self.assertEqual(actual, """IS_name\tsample_id\tcontig\titr1_start_position\titr1_end_position\titr2_start_position\titr2_end_position\titr_cluster\nIS_cluster_0_length_1021\ttest\tNODE_823_length_1805_cov_1014.02\t246\t271\t1242\t1266\t0\n""")

        fasta_name = self.TEST_OUTPUT_PREFIX + '_candidate_insertion_sequences.fasta'
        fasta = open(fasta_name, "r")
        actual = "".join(fasta.readlines())
        os.remove(fasta_name)
        self.assertEqual(actual, """>IS_cluster_0_length_1021\nAATATTGTTTTACATCCTGCATCTTACATTTATATGTATTATGATATATAACAAAGGTACTATATCATATATGATATATAATATATATATTTGTAGTACATATTATAATTTTTATATATTATAATTTATATAATGATAAAATTTTTTATAGAATATAAATAATTATATATGATTATATAATTCTACATATTATAAATGAAAATATGTAGAATTACATGTTTATCTAGTTATATGTATAAAAATATAATTATATATATAGTTATATATATCACATAATAGACAACATACATATTATATAATGTGTAATATATATTATAGAATATATGATACATGATATATATGATATTTAATATAACATAATATAATCTATATTATAATATTATGTATGTTACTTATATTGGGTGATATGTAATATATATTATGTAATATGAAATAATATAATATATATTATATTATGATATTTTATGTAAATTATGTTATAAAAGTATATATAACATAATATATAGTTATATATAATATATTATATAGTTATATATATTATTATATTATAATTTATGACATATATAACAGTTATATATAATTTCAAATAGTTATATATAACAGAATATAAAACATATAGAATACATATCATAAAATATATATTGTATACCATATATATTAGGTATCATATATTGTATGTTATATATCATATACTGTATATCATATATCATATATTATATGTTATATATATCACATTCCACATATGGTATATTATAAATTATACATAATATATTGCATTCTTTATTTTACATGTAATAAATTATACATTATAGGTAACATAGTATATATTGTATGTAACCGGTATATTCTATGTAAAATATATAATATATAACACATGTGACATGTAATTAATAACATAGTCATATAATATATATTTAATATAGGTATTATACATACAACGTATTATATATAATATATAATATATATTACATATGTTATAAAATATAATATATATTATAAGATGCAGGATGTAAAACAATAT\n""")

        for blast_file in Path('./').glob('blastn_*'):
            os.remove(blast_file)

        for itr_fasta in Path('./').glob('itr*_tmp.fasta'):
            os.remove(itr_fasta)

    def test_arguments(self):
        actual = get_arguments().parse_args(
            ['--cdhit_cluster_file', 'clstr_file', '--info_tab_file', 'tab_file',
            '--assemblies_fasta_file', 'assembly_file', '--min_is_len', '500',
            '--max_is_len', '3000', '--min_itr_len', '25', '--max_itr_len', '50', '--cpus', '8', '--output_prefix', 'output_prefix'])
        self.assertEqual(actual, argparse.Namespace(cluster_file='clstr_file', tab_file='tab_file', assemblies_file='assembly_file', min_is_len = 500, max_is_len = 3000, min_itr_len = 25, max_itr_len = 50, cpus = 8, output_prefix='output_prefix'))

    @patch('bin.assign_ITRs.create_assembly_bins')
    @patch('bin.assign_ITRs.create_cluster_dictionary')
    @patch('bin.assign_ITRs.bin_positions')
    @patch('bin.assign_ITRs.write_itr_annotations')
    def test_main(self, mock_write_itr_annotations, mock_bin_positions, mock_create_cluster_dictionary, mock_create_assembly_bins):
        args = get_arguments().parse_args(
            ['--cdhit_cluster_file', 'clstr_file', '--info_tab_file', 'tab_file',
            '--assemblies_fasta_file', 'assembly_file', '--min_is_len', '500', '--max_is_len', '3000', '--min_itr_len', '25', '--max_itr_len', '50', '--cpus', '8', '--output_prefix', 'output_prefix'])
        mock_create_assembly_bins.return_value = 'assembly_bins_dict'
        mock_create_cluster_dictionary.return_value = 'cl_dict'
        mock_bin_positions.return_value = 'clusters_positions'

        main(args)

        self.assertEqual(mock_create_assembly_bins.call_args_list, [call('assembly_file')])
        self.assertEqual(mock_create_cluster_dictionary.call_args_list, [call('clstr_file')])
        self.assertEqual(mock_bin_positions.call_args_list, [call('cl_dict', 'tab_file', 'assembly_bins_dict', 'output_prefix')])
        self.assertEqual(mock_write_itr_annotations.call_args_list, [call('clusters_positions', 'assembly_bins_dict', 500, 3000, 25, 50, 'output_prefix', 8)])
