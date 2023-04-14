import argparse, os
import unittest
from unittest.mock import patch, call, ANY

from bin.get_insertion_sequences import *

class TestGetInsertionSequences(unittest.TestCase):
    TEST_INTERPROSCAN_FILE = 'tests/data/input/test_scaffolds_filtered_transposase.tsv'
    TEST_CONTIG_INFO = 'tests/data/input/test_candidate_insertion_sequences_info.txt'
    TEST_CONTIG_FASTA = 'tests/data/input/test2_candidate_insertion_sequences.fasta'
    TEST_OUTPUT_PREFIX = 'tests/data/output/test'
    TEST_OUTPUT_PREFIX_MAIN = 'tests/data/output/test2'
    TEST_OUTPUT_PREFIX_MAIN2 = 'tests/data/output/test3'

    def test_get_contig_info(self):

        contig_transposase_info = contigTransposaseInfo()
        contig_transposase_info.load_contig_info(self.TEST_INTERPROSCAN_FILE)

        self.maxDiff = None
        self.assertEqual(contig_transposase_info.get_contig_info(), {'NODE_1_length_49561_cov_5.4235_1': [('IPR036515',
                                       'Transposase IS200-like superfamily',
                                       49,
                                       549),
                                      ('IPR036515',
                                       'Transposase IS200-like superfamily',
                                       52,
                                       534)],
 'NODE_1_length_49561_cov_5.4235_2': [('PTHR36966',
                                       'REP-ASSOCIATED TYROSINE TRANSPOSASE',
                                       1,
                                       528),
                                      ('IPR002686',
                                       'Transposase IS200-like',
                                       67,
                                       522),
                                      ('IPR002686',
                                       'Transposase IS200-like',
                                       82,
                                       513)]}
        )

    def test_write_insertion_sequence_info(self):

        contig_transposase_info = contigTransposaseInfo()
        contig_transposase_info.load_contig_info(self.TEST_INTERPROSCAN_FILE)

        actual = write_insertion_sequences_info(contig_transposase_info, self.TEST_CONTIG_INFO, self.TEST_OUTPUT_PREFIX)

        self.maxDiff = None
        self.assertEqual(actual, {'NODE_1_length_49561_cov_5.4235_1': 'IS_tests/data/output/test_NODE_1_length_49561_cov_5.4235_1_10_570',
'NODE_1_length_49561_cov_5.4235_2': 'IS_tests/data/output/test_NODE_1_length_49561_cov_5.4235_2_10_600'}
        )

        txt = open(f'{self.TEST_OUTPUT_PREFIX}_insertion_sequences_info.txt', "r")
        content = "".join(txt.readlines())
        os.remove(f'{self.TEST_OUTPUT_PREFIX}_insertion_sequences_info.txt')
        self.assertEqual(content,"""IS_name	sample_id	contig	is_start	is_end	description
IS_tests/data/output/test_NODE_1_length_49561_cov_5.4235_1_10_570	tests/data/output/test	NODE_1_length_49561_cov_5.4235_1	10	570	IPR036515:Transposase IS200-like superfamily:49:549;IPR036515:Transposase IS200-like superfamily:52:534
IS_tests/data/output/test_NODE_1_length_49561_cov_5.4235_2_10_600	tests/data/output/test	NODE_1_length_49561_cov_5.4235_2	10	600	PTHR36966:REP-ASSOCIATED TYROSINE TRANSPOSASE:1:528;IPR002686:Transposase IS200-like:67:522;IPR002686:Transposase IS200-like:82:513
""")

    def test_write_insertion_sequence_fasta(self):

        contig_is_info = {'NODE_1_length_49561_cov_5.4235_1': 'IS_test_NODE_1_length_49561_cov_5.4235_1_10_570'}

        write_insertion_sequences_fasta(contig_is_info, self.TEST_CONTIG_FASTA, self.TEST_OUTPUT_PREFIX)

        fasta = open(f'{self.TEST_OUTPUT_PREFIX}_insertion_sequences.fasta', "r")
        content = "".join(fasta.readlines())
        os.remove(f'{self.TEST_OUTPUT_PREFIX}_insertion_sequences.fasta')
        self.assertEqual(content,""">IS_test_NODE_1_length_49561_cov_5.4235_1_10_570
GTACGGTCTAGGAAGCACCTTTACCCTTGTTCTCAATCTCTCTGGTAGTGAAAATAAAGCTTGTAAGTGTCTAGTGTTTTGAGGTAAATAGAAAATATTGAGCTAACGGAAAATTCCGTTAGTTCTTTTTTGTATTTTAGGAGAAAAAATGTCAGATAATATTGAAATTGTTATGGATATTGATTCGTTAGAAATTTTCGATGTAATAGATCTTGATAAGAATGAAAATATGGTTGATTGTATCCCTGAAGGAAATTCTGAAGAAAAGGTATCTCCATACGATGACCTCAATAGGGTTTGGCTTAGAAACGGGTATTATTATAAGCCGTCTAGTCTACTTTGGGAAGTTATTAAAGATAATGTTATTCATGAATTGTATCCTAGAAATAAACATGGCAAATATCCTTACTACAGAAGAGCAAGAGACGGAACATGGGGACAATCTAGCCAATTGACTCAGGATGAAGACTTGTGGTTGACTCTAGAGTTTAATCCAGAATTGATGGGGTATGTGTAGGAAGGAGTAAAGAATGGAAACTGGAATTTGTATACGATGTGCTT
""")

    def test_arguments(self):
        actual = get_arguments().parse_args(
            ['--contig_fasta', 'contig_fasta', '--contig_info', 'contig_info', '--interproscan_info', 'interproscan_info', '--output_prefix', 'output_prefix'])
        self.assertEqual(actual, argparse.Namespace(contig_fasta='contig_fasta', contig_info='contig_info', interpro_info='interproscan_info', output_prefix='output_prefix'))

    def test_main(self):
        args = get_arguments().parse_args(
            ['--contig_fasta', self.TEST_CONTIG_FASTA, '--contig_info', self.TEST_CONTIG_INFO, '--interproscan_info', self.TEST_INTERPROSCAN_FILE, '--output_prefix', self.TEST_OUTPUT_PREFIX_MAIN])

        main(args)

        txt = open(f'{self.TEST_OUTPUT_PREFIX_MAIN}_insertion_sequences_info.txt', "r")
        content = "".join(txt.readlines())
        os.remove(f'{self.TEST_OUTPUT_PREFIX_MAIN}_insertion_sequences_info.txt')
        self.assertEqual(content,"""IS_name	sample_id	contig	is_start	is_end	description
IS_tests/data/output/test2_NODE_1_length_49561_cov_5.4235_1_10_570	tests/data/output/test2	NODE_1_length_49561_cov_5.4235_1	10	570	IPR036515:Transposase IS200-like superfamily:49:549;IPR036515:Transposase IS200-like superfamily:52:534
IS_tests/data/output/test2_NODE_1_length_49561_cov_5.4235_2_10_600	tests/data/output/test2	NODE_1_length_49561_cov_5.4235_2	10	600	PTHR36966:REP-ASSOCIATED TYROSINE TRANSPOSASE:1:528;IPR002686:Transposase IS200-like:67:522;IPR002686:Transposase IS200-like:82:513
""")

        fasta = open(f'{self.TEST_OUTPUT_PREFIX_MAIN}_insertion_sequences.fasta', "r")
        content = "".join(fasta.readlines())
        self.maxDiff = None
        os.remove(f'{self.TEST_OUTPUT_PREFIX_MAIN}_insertion_sequences.fasta')
        self.assertEqual(content,""">IS_tests/data/output/test2_NODE_1_length_49561_cov_5.4235_1_10_570
GTACGGTCTAGGAAGCACCTTTACCCTTGTTCTCAATCTCTCTGGTAGTGAAAATAAAGCTTGTAAGTGTCTAGTGTTTTGAGGTAAATAGAAAATATTGAGCTAACGGAAAATTCCGTTAGTTCTTTTTTGTATTTTAGGAGAAAAAATGTCAGATAATATTGAAATTGTTATGGATATTGATTCGTTAGAAATTTTCGATGTAATAGATCTTGATAAGAATGAAAATATGGTTGATTGTATCCCTGAAGGAAATTCTGAAGAAAAGGTATCTCCATACGATGACCTCAATAGGGTTTGGCTTAGAAACGGGTATTATTATAAGCCGTCTAGTCTACTTTGGGAAGTTATTAAAGATAATGTTATTCATGAATTGTATCCTAGAAATAAACATGGCAAATATCCTTACTACAGAAGAGCAAGAGACGGAACATGGGGACAATCTAGCCAATTGACTCAGGATGAAGACTTGTGGTTGACTCTAGAGTTTAATCCAGAATTGATGGGGTATGTGTAGGAAGGAGTAAAGAATGGAAACTGGAATTTGTATACGATGTGCTT
>IS_tests/data/output/test2_NODE_1_length_49561_cov_5.4235_2_10_600
GTACGGTCTAGGAAGCACCTTTACCCTTGTTCTCAATCTCTCTGGTAGTGAAAATAAAGCTTGTAAGTGTCTAGTGTTTTGAGGTAAATAGAAAATATTGAGCTAACGGAAAATTCCGTTAGTTCTTTTTTGTATTTTAGGAGAAAAAATGTCAGATAATATTGAAATTGTTATGGATATTGATTCGTTAGAAATTTTCGATGTAATAGATCTTGATAAGAATGAAAATATGGTTGATTGTATCCCTGAAGGAAATTCTGAAGAAAAGGTATCTCCATACGATGACCTCAATAGGGTTTGGCTTAGAAACGGGTATTATTATAAGCCGTCTAGTCTACTTTGGGAAGTTATTAAAGATAATGTTATTCATGAATTGTATCCTAGAAATAAACATGGCAAATATCCTTACTACAGAAGAGCAAGAGACGGAACATGGGGACAATCTAGCCAATTGACTCAGGATGAAGACTTGTGGTTGACTCTAGAGTTTAATCCAGAATTGATGGGGTATGTGTAGGAAGGAGTAAAGAATGGAAACTGGAATTTGTATACGATGTGCTTTCAACTGGGTTACACCTTGCATCAATGAAA
""")
