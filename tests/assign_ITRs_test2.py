import argparse, os
import unittest
from unittest.mock import patch, call, ANY

from bin.assign_ITRs import *

class TestAssignITRs2(unittest.TestCase):
    TEST_CLSTR = 'tests/data/input/test4.fasta.clstr'
    TEST_INFO_TAB = 'tests/data/input/test4_contigs_reads_ir_position_info.tab'
    TEST_OUTPUT_PREFIX = 'tests/data/output/test3'
    TEST_CONTIG_FASTA = 'tests/data/input/test3_assemblies_filtered.fasta'

    def test_assembly_bins(self):
        actual = create_assembly_bins(self.TEST_CONTIG_FASTA)

        self.assertEqual(len(actual['ERR1624733.20280_5_1.13']), 26459)

    def test_create_cluster_dictionary(self):
        actual = create_cluster_dictionary(self.TEST_CLSTR)
        self.assertEqual(actual[1][95], '')
        self.assertEqual(actual[0][95], '')
        self.assertEqual(actual[0][228004], '1116')

    def test_bin_positions(self):
        assembly_bins_dict = create_assembly_bins(self.TEST_CONTIG_FASTA)
        cl_dict = create_cluster_dictionary(self.TEST_CLSTR)
        actual = bin_positions(cl_dict, self.TEST_INFO_TAB, assembly_bins_dict, self.TEST_OUTPUT_PREFIX)
        self.assertEqual(actual['ERR1624733.20280_5_1.13']['1116'][0], '0')
        self.assertEqual(actual['ERR1624733.20280_5_1.13']['1116'][11653], '1')

    def test_count_bins(self):
        assembly_bins_dict = create_assembly_bins(self.TEST_CONTIG_FASTA)
        cl_dict = create_cluster_dictionary(self.TEST_CLSTR)
        bins = bin_positions(cl_dict, self.TEST_INFO_TAB, assembly_bins_dict, self.TEST_OUTPUT_PREFIX)
        actual = count_bins(bins['ERR1624733.20280_5_1.13']['1116'])

        self.assertEqual(actual, [('0', 11630), ('1', 77), ('0', 2594), ('1', 89), ('0', 12069)])

    def test_get_itrs_from_count_bins_101(self):
        assembly_bins_dict = create_assembly_bins(self.TEST_CONTIG_FASTA)
        cl_dict = create_cluster_dictionary(self.TEST_CLSTR)
        clusters_positions = bin_positions(cl_dict, self.TEST_INFO_TAB, assembly_bins_dict, self.TEST_OUTPUT_PREFIX)
        count_bins_out = count_bins(clusters_positions['ERR1624733.20280_5_1.13']['1116'])
        actual = get_itrs_from_count_bins(count_bins_out, 500, 3000, 25, 50)

        self.assertEqual(actual, [])
