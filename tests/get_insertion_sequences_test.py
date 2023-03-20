import argparse, os
import unittest
from unittest.mock import patch, call, ANY

from bin.get_insertion_sequences import *

class TestCandidateITRReadsAndISContigs(unittest.TestCase):
    TEST_CLSTR = 'tests/data/input/test.fasta.clstr'
    TEST_CONTIG_FILE = 'tests/data/input/test_scaffolds_filtered_transposase.tsv'
    TEST_SAM_FILE = 'tests/data/input/test.sam.mapped.sorted'
    TEST_CONTIG_FASTA = 'tests/data/input/test_assemblies_filtered.fasta'
    TEST_OUTPUT_PREFIX = 'tests/data/output/test'
    TEST_OUTPUT_PREFIX_MAIN = 'tests/data/output/test2'
    TEST_OUTPUT_PREFIX_MAIN2 = 'tests/data/output/test3'

    def test_get_ir_cluster(self):

        clusterInfo = irClusterInfo()
        clusterInfo.load_ir_clusters(self.TEST_CLSTR)

        self.assertEqual(clusterInfo.get_ir_cluster("Seq7284_nstart_SRS011144_nend_SRR062032.384055.1_f1_LCoord_34_RCoord_60"), '0')
        self.assertEqual(clusterInfo.get_ir_cluster("Seq632295_nstart_SRS011144_nend_SRR062032.27504454.1_f1_LCoord_47_RCoord_72"), '8601')

    def test_get_contig_info(self):

        contig_transposase_info = contigTransposaseInfo()
        contig_transposase_info.load_contig_info(self.TEST_CONTIG_FILE)

        self.maxDiff = None
        self.assertEqual(contig_transposase_info.get_contig_info(), {'NODE_1_length_49561_cov_5.4235': [['IPR036515',
                                      'Transposase IS200-like superfamily',
                                      49,
                                      549],
                                     ['IPR036515',
                                      'Transposase IS200-like superfamily',
                                      52,
                                      534],
                                     ['PTHR36966',
                                      'REP-ASSOCIATED TYROSINE TRANSPOSASE',
                                      1,
                                      528],
                                     ['IPR002686',
                                      'Transposase IS200-like',
                                      67,
                                      522],
                                     ['IPR002686',
                                      'Transposase IS200-like',
                                      82,
                                      513]],
  'NODE_5_length_32495_cov_4.56202': [['PTHR35004',
                                       'TRANSPOSASE RV3428C-RELATED',
                                       1,
                                       864],
                                      ['IPR012337',
                                       'Ribonuclease H-like superfamily',
                                       346,
                                       837],
                                      ['IPR036397',
                                       'Ribonuclease H superfamily',
                                       346,
                                       837]],
  'NODE_6_length_32110_cov_3.8447': [['IPR013762',
                                      'Integrase-like, catalytic domain '
                                      'superfamily',
                                      448,
                                      975]]})

    def test_check_irs_flanking(self):

        contig_transposase_info = contigTransposaseInfo()
        contig_transposase_info.load_contig_info(self.TEST_CONTIG_FILE)
        actual = contig_transposase_info.check_irs_flanking('NODE_1_length_49561_cov_5.4235', [[10, 30], [550,570], [0,20], [580,600]], 300, 5000)
        self.maxDiff = None
        self.assertEqual(actual, {(0, 570): [['IPR036515', 'Transposase IS200-like superfamily', 49, 549],
            ['IPR036515', 'Transposase IS200-like superfamily', 52, 534],
            ['IPR002686', 'Transposase IS200-like', 67, 522],
            ['IPR002686', 'Transposase IS200-like', 82, 513]],
 (0, 600): [['IPR036515', 'Transposase IS200-like superfamily', 49, 549],
            ['IPR036515', 'Transposase IS200-like superfamily', 52, 534],
            ['IPR002686', 'Transposase IS200-like', 67, 522],
            ['IPR002686', 'Transposase IS200-like', 82, 513]],
 (10, 570): [['IPR036515', 'Transposase IS200-like superfamily', 49, 549],
             ['IPR036515', 'Transposase IS200-like superfamily', 52, 534],
             ['IPR002686', 'Transposase IS200-like', 67, 522],
             ['IPR002686', 'Transposase IS200-like', 82, 513]],
 (10, 600): [['IPR036515', 'Transposase IS200-like superfamily', 49, 549],
             ['IPR036515', 'Transposase IS200-like superfamily', 52, 534],
             ['IPR002686', 'Transposase IS200-like', 67, 522],
             ['IPR002686', 'Transposase IS200-like', 82, 513]]}
        )

    def test_process_sam_file(self):

        clusterInfo = irClusterInfo()
        clusterInfo.load_ir_clusters(self.TEST_CLSTR)

        actual = process_sam_file(self.TEST_SAM_FILE, clusterInfo)

        self.assertEqual(actual, {
'NODE_1_length_49561_cov_5.4235': {'0': [[48, 147], [546, 645], [531, 630]]},
'NODE_5_length_32495_cov_4.56202': {'8601': [[51, 150]]}})

    def test_write_insertion_sequence_info(self):

        clusterInfo = irClusterInfo()
        clusterInfo.load_ir_clusters(self.TEST_CLSTR)
        mapping_info = {
'NODE_1_length_49561_cov_5.4235': {'0': [[10, 30], [550,570], [1,20], [580,600]]},
'NODE_5_length_32495_cov_4.56202': {'8601': [[51, 150]]}}
        contig_transposase_info = contigTransposaseInfo()
        contig_transposase_info.load_contig_info(self.TEST_CONTIG_FILE)

        actual = write_insertion_sequences_info(mapping_info, contig_transposase_info, self.TEST_OUTPUT_PREFIX, 500, 3000)

        self.maxDiff = None
        self.assertEqual(actual, {'NODE_1_length_49561_cov_5.4235': [['IS_length_561_IPR036515_49_549_IPR036515_52_534_IPR002686_67_522_IPR002686_82_513',
                                      10,
                                      570],
                                     ['IS_length_591_IPR036515_49_549_IPR036515_52_534_IPR002686_67_522_IPR002686_82_513',
                                      10,
                                      600],
                                     ['IS_length_570_IPR036515_49_549_IPR036515_52_534_IPR002686_67_522_IPR002686_82_513',
                                      1,
                                      570],
                                     ['IS_length_600_IPR036515_49_549_IPR036515_52_534_IPR002686_67_522_IPR002686_82_513',
                                      1,
                                      600]]})

        txt = open(f'{self.TEST_OUTPUT_PREFIX}_insertion_sequences_info.txt', "r")
        content = "".join(txt.readlines())
        os.remove(f'{self.TEST_OUTPUT_PREFIX}_insertion_sequences_info.txt')
        self.assertEqual(content,"""IS_name	sample_id	contig	is_start	is_end	description
IS_length_561_IPR036515_49_549_IPR036515_52_534_IPR002686_67_522_IPR002686_82_513	tests/data/output/test	NODE_1_length_49561_cov_5.4235	10	570	Transposase IS200-like superfamily;Transposase IS200-like superfamily;Transposase IS200-like;Transposase IS200-like
IS_length_591_IPR036515_49_549_IPR036515_52_534_IPR002686_67_522_IPR002686_82_513	tests/data/output/test	NODE_1_length_49561_cov_5.4235	10	600	Transposase IS200-like superfamily;Transposase IS200-like superfamily;Transposase IS200-like;Transposase IS200-like
IS_length_570_IPR036515_49_549_IPR036515_52_534_IPR002686_67_522_IPR002686_82_513	tests/data/output/test	NODE_1_length_49561_cov_5.4235	1	570	Transposase IS200-like superfamily;Transposase IS200-like superfamily;Transposase IS200-like;Transposase IS200-like
IS_length_600_IPR036515_49_549_IPR036515_52_534_IPR002686_67_522_IPR002686_82_513	tests/data/output/test	NODE_1_length_49561_cov_5.4235	1	600	Transposase IS200-like superfamily;Transposase IS200-like superfamily;Transposase IS200-like;Transposase IS200-like
""")

    def test_write_insertion_sequence_fasta(self):

        contig_is_info = {'NODE_1_length_49561_cov_5.4235': [['IS_length_561_IPR036515_49_549_IPR036515_52_534_IPR002686_67_522_IPR002686_82_513',
                                      10,
                                      570],
                                     ['IS_length_591_IPR036515_49_549_IPR036515_52_534_IPR002686_67_522_IPR002686_82_513',
                                      10,
                                      600],
                                     ['IS_length_570_IPR036515_49_549_IPR036515_52_534_IPR002686_67_522_IPR002686_82_513',
                                      1,
                                      570],
                                     ['IS_length_600_IPR036515_49_549_IPR036515_52_534_IPR002686_67_522_IPR002686_82_513',
                                      1,
                                      600]]}

        write_insertion_sequences_fasta(contig_is_info, self.TEST_CONTIG_FASTA, self.TEST_OUTPUT_PREFIX)

        fasta = open(f'{self.TEST_OUTPUT_PREFIX}_insertion_sequences.fasta', "r")
        content = "".join(fasta.readlines())
        os.remove(f'{self.TEST_OUTPUT_PREFIX}_insertion_sequences.fasta')
        self.assertEqual(content,""">IS_length_561_IPR036515_49_549_IPR036515_52_534_IPR002686_67_522_IPR002686_82_513
GTACGGTCTAGGAAGCACCTTTACCCTTGTTCTCAATCTCTCTGGTAGTGAAAATAAAGCTTGTAAGTGTCTAGTGTTTTGAGGTAAATAGAAAATATTGAGCTAACGGAAAATTCCGTTAGTTCTTTTTTGTATTTTAGGAGAAAAAATGTCAGATAATATTGAAATTGTTATGGATATTGATTCGTTAGAAATTTTCGATGTAATAGATCTTGATAAGAATGAAAATATGGTTGATTGTATCCCTGAAGGAAATTCTGAAGAAAAGGTATCTCCATACGATGACCTCAATAGGGTTTGGCTTAGAAACGGGTATTATTATAAGCCGTCTAGTCTACTTTGGGAAGTTATTAAAGATAATGTTATTCATGAATTGTATCCTAGAAATAAACATGGCAAATATCCTTACTACAGAAGAGCAAGAGACGGAACATGGGGACAATCTAGCCAATTGACTCAGGATGAAGACTTGTGGTTGACTCTAGAGTTTAATCCAGAATTGATGGGGTATGTGTAGGAAGGAGTAAAGAATGGAAACTGGAATTTGTATACGATGTGCTT
>IS_length_591_IPR036515_49_549_IPR036515_52_534_IPR002686_67_522_IPR002686_82_513
GTACGGTCTAGGAAGCACCTTTACCCTTGTTCTCAATCTCTCTGGTAGTGAAAATAAAGCTTGTAAGTGTCTAGTGTTTTGAGGTAAATAGAAAATATTGAGCTAACGGAAAATTCCGTTAGTTCTTTTTTGTATTTTAGGAGAAAAAATGTCAGATAATATTGAAATTGTTATGGATATTGATTCGTTAGAAATTTTCGATGTAATAGATCTTGATAAGAATGAAAATATGGTTGATTGTATCCCTGAAGGAAATTCTGAAGAAAAGGTATCTCCATACGATGACCTCAATAGGGTTTGGCTTAGAAACGGGTATTATTATAAGCCGTCTAGTCTACTTTGGGAAGTTATTAAAGATAATGTTATTCATGAATTGTATCCTAGAAATAAACATGGCAAATATCCTTACTACAGAAGAGCAAGAGACGGAACATGGGGACAATCTAGCCAATTGACTCAGGATGAAGACTTGTGGTTGACTCTAGAGTTTAATCCAGAATTGATGGGGTATGTGTAGGAAGGAGTAAAGAATGGAAACTGGAATTTGTATACGATGTGCTTTCAACTGGGTTACACCTTGCATCAATGAAA
>IS_length_570_IPR036515_49_549_IPR036515_52_534_IPR002686_67_522_IPR002686_82_513
CAGCAGCCAGTACGGTCTAGGAAGCACCTTTACCCTTGTTCTCAATCTCTCTGGTAGTGAAAATAAAGCTTGTAAGTGTCTAGTGTTTTGAGGTAAATAGAAAATATTGAGCTAACGGAAAATTCCGTTAGTTCTTTTTTGTATTTTAGGAGAAAAAATGTCAGATAATATTGAAATTGTTATGGATATTGATTCGTTAGAAATTTTCGATGTAATAGATCTTGATAAGAATGAAAATATGGTTGATTGTATCCCTGAAGGAAATTCTGAAGAAAAGGTATCTCCATACGATGACCTCAATAGGGTTTGGCTTAGAAACGGGTATTATTATAAGCCGTCTAGTCTACTTTGGGAAGTTATTAAAGATAATGTTATTCATGAATTGTATCCTAGAAATAAACATGGCAAATATCCTTACTACAGAAGAGCAAGAGACGGAACATGGGGACAATCTAGCCAATTGACTCAGGATGAAGACTTGTGGTTGACTCTAGAGTTTAATCCAGAATTGATGGGGTATGTGTAGGAAGGAGTAAAGAATGGAAACTGGAATTTGTATACGATGTGCTT
>IS_length_600_IPR036515_49_549_IPR036515_52_534_IPR002686_67_522_IPR002686_82_513
CAGCAGCCAGTACGGTCTAGGAAGCACCTTTACCCTTGTTCTCAATCTCTCTGGTAGTGAAAATAAAGCTTGTAAGTGTCTAGTGTTTTGAGGTAAATAGAAAATATTGAGCTAACGGAAAATTCCGTTAGTTCTTTTTTGTATTTTAGGAGAAAAAATGTCAGATAATATTGAAATTGTTATGGATATTGATTCGTTAGAAATTTTCGATGTAATAGATCTTGATAAGAATGAAAATATGGTTGATTGTATCCCTGAAGGAAATTCTGAAGAAAAGGTATCTCCATACGATGACCTCAATAGGGTTTGGCTTAGAAACGGGTATTATTATAAGCCGTCTAGTCTACTTTGGGAAGTTATTAAAGATAATGTTATTCATGAATTGTATCCTAGAAATAAACATGGCAAATATCCTTACTACAGAAGAGCAAGAGACGGAACATGGGGACAATCTAGCCAATTGACTCAGGATGAAGACTTGTGGTTGACTCTAGAGTTTAATCCAGAATTGATGGGGTATGTGTAGGAAGGAGTAAAGAATGGAAACTGGAATTTGTATACGATGTGCTTTCAACTGGGTTACACCTTGCATCAATGAAA
""")

    def test_arguments(self):
        actual = get_arguments().parse_args(
            ['--contig_fasta', 'contig_fasta', '--contig_info', 'contig_info', '--sam_file', 'sam_file', '--clstr_file', 'clstr_file', '--min_is_len', '500', '--max_is_len', '3000', '--output_prefix', 'output_prefix'])
        self.assertEqual(actual, argparse.Namespace(contig_fasta='contig_fasta', contig_info='contig_info', sam_file='sam_file', clstr_file='clstr_file', min_is_len=500, max_is_len=3000, output_prefix='output_prefix'))

    @patch('bin.get_insertion_sequences.process_sam_file')
    def test_main(self, mock_process_sam_file):
        args = get_arguments().parse_args(
            ['--contig_fasta', self.TEST_CONTIG_FASTA, '--contig_info', self.TEST_CONTIG_FILE, '--sam_file', 'sam_file', '--clstr_file', self.TEST_CLSTR, '--output_prefix', self.TEST_OUTPUT_PREFIX_MAIN])

        mock_process_sam_file.return_value = mapping_info = {
'NODE_1_length_49561_cov_5.4235': {'0': [[10, 30], [550,570], [1,20], [580,600]]},
'NODE_5_length_32495_cov_4.56202': {'8601': [[51, 150]]}}

        main(args)

        txt = open(f'{self.TEST_OUTPUT_PREFIX_MAIN}_insertion_sequences_info.txt', "r")
        content = "".join(txt.readlines())
        os.remove(f'{self.TEST_OUTPUT_PREFIX_MAIN}_insertion_sequences_info.txt')
        self.assertEqual(content,"""IS_name	sample_id	contig	is_start	is_end	description
IS_length_561_IPR036515_49_549_IPR036515_52_534_IPR002686_67_522_IPR002686_82_513	tests/data/output/test2	NODE_1_length_49561_cov_5.4235	10	570	Transposase IS200-like superfamily;Transposase IS200-like superfamily;Transposase IS200-like;Transposase IS200-like
IS_length_591_IPR036515_49_549_IPR036515_52_534_IPR002686_67_522_IPR002686_82_513	tests/data/output/test2	NODE_1_length_49561_cov_5.4235	10	600	Transposase IS200-like superfamily;Transposase IS200-like superfamily;Transposase IS200-like;Transposase IS200-like
IS_length_570_IPR036515_49_549_IPR036515_52_534_IPR002686_67_522_IPR002686_82_513	tests/data/output/test2	NODE_1_length_49561_cov_5.4235	1	570	Transposase IS200-like superfamily;Transposase IS200-like superfamily;Transposase IS200-like;Transposase IS200-like
IS_length_600_IPR036515_49_549_IPR036515_52_534_IPR002686_67_522_IPR002686_82_513	tests/data/output/test2	NODE_1_length_49561_cov_5.4235	1	600	Transposase IS200-like superfamily;Transposase IS200-like superfamily;Transposase IS200-like;Transposase IS200-like
""")

        fasta = open(f'{self.TEST_OUTPUT_PREFIX_MAIN}_insertion_sequences.fasta', "r")
        content = "".join(fasta.readlines())
        os.remove(f'{self.TEST_OUTPUT_PREFIX_MAIN}_insertion_sequences.fasta')
        self.assertEqual(content,""">IS_length_561_IPR036515_49_549_IPR036515_52_534_IPR002686_67_522_IPR002686_82_513
GTACGGTCTAGGAAGCACCTTTACCCTTGTTCTCAATCTCTCTGGTAGTGAAAATAAAGCTTGTAAGTGTCTAGTGTTTTGAGGTAAATAGAAAATATTGAGCTAACGGAAAATTCCGTTAGTTCTTTTTTGTATTTTAGGAGAAAAAATGTCAGATAATATTGAAATTGTTATGGATATTGATTCGTTAGAAATTTTCGATGTAATAGATCTTGATAAGAATGAAAATATGGTTGATTGTATCCCTGAAGGAAATTCTGAAGAAAAGGTATCTCCATACGATGACCTCAATAGGGTTTGGCTTAGAAACGGGTATTATTATAAGCCGTCTAGTCTACTTTGGGAAGTTATTAAAGATAATGTTATTCATGAATTGTATCCTAGAAATAAACATGGCAAATATCCTTACTACAGAAGAGCAAGAGACGGAACATGGGGACAATCTAGCCAATTGACTCAGGATGAAGACTTGTGGTTGACTCTAGAGTTTAATCCAGAATTGATGGGGTATGTGTAGGAAGGAGTAAAGAATGGAAACTGGAATTTGTATACGATGTGCTT
>IS_length_591_IPR036515_49_549_IPR036515_52_534_IPR002686_67_522_IPR002686_82_513
GTACGGTCTAGGAAGCACCTTTACCCTTGTTCTCAATCTCTCTGGTAGTGAAAATAAAGCTTGTAAGTGTCTAGTGTTTTGAGGTAAATAGAAAATATTGAGCTAACGGAAAATTCCGTTAGTTCTTTTTTGTATTTTAGGAGAAAAAATGTCAGATAATATTGAAATTGTTATGGATATTGATTCGTTAGAAATTTTCGATGTAATAGATCTTGATAAGAATGAAAATATGGTTGATTGTATCCCTGAAGGAAATTCTGAAGAAAAGGTATCTCCATACGATGACCTCAATAGGGTTTGGCTTAGAAACGGGTATTATTATAAGCCGTCTAGTCTACTTTGGGAAGTTATTAAAGATAATGTTATTCATGAATTGTATCCTAGAAATAAACATGGCAAATATCCTTACTACAGAAGAGCAAGAGACGGAACATGGGGACAATCTAGCCAATTGACTCAGGATGAAGACTTGTGGTTGACTCTAGAGTTTAATCCAGAATTGATGGGGTATGTGTAGGAAGGAGTAAAGAATGGAAACTGGAATTTGTATACGATGTGCTTTCAACTGGGTTACACCTTGCATCAATGAAA
>IS_length_570_IPR036515_49_549_IPR036515_52_534_IPR002686_67_522_IPR002686_82_513
CAGCAGCCAGTACGGTCTAGGAAGCACCTTTACCCTTGTTCTCAATCTCTCTGGTAGTGAAAATAAAGCTTGTAAGTGTCTAGTGTTTTGAGGTAAATAGAAAATATTGAGCTAACGGAAAATTCCGTTAGTTCTTTTTTGTATTTTAGGAGAAAAAATGTCAGATAATATTGAAATTGTTATGGATATTGATTCGTTAGAAATTTTCGATGTAATAGATCTTGATAAGAATGAAAATATGGTTGATTGTATCCCTGAAGGAAATTCTGAAGAAAAGGTATCTCCATACGATGACCTCAATAGGGTTTGGCTTAGAAACGGGTATTATTATAAGCCGTCTAGTCTACTTTGGGAAGTTATTAAAGATAATGTTATTCATGAATTGTATCCTAGAAATAAACATGGCAAATATCCTTACTACAGAAGAGCAAGAGACGGAACATGGGGACAATCTAGCCAATTGACTCAGGATGAAGACTTGTGGTTGACTCTAGAGTTTAATCCAGAATTGATGGGGTATGTGTAGGAAGGAGTAAAGAATGGAAACTGGAATTTGTATACGATGTGCTT
>IS_length_600_IPR036515_49_549_IPR036515_52_534_IPR002686_67_522_IPR002686_82_513
CAGCAGCCAGTACGGTCTAGGAAGCACCTTTACCCTTGTTCTCAATCTCTCTGGTAGTGAAAATAAAGCTTGTAAGTGTCTAGTGTTTTGAGGTAAATAGAAAATATTGAGCTAACGGAAAATTCCGTTAGTTCTTTTTTGTATTTTAGGAGAAAAAATGTCAGATAATATTGAAATTGTTATGGATATTGATTCGTTAGAAATTTTCGATGTAATAGATCTTGATAAGAATGAAAATATGGTTGATTGTATCCCTGAAGGAAATTCTGAAGAAAAGGTATCTCCATACGATGACCTCAATAGGGTTTGGCTTAGAAACGGGTATTATTATAAGCCGTCTAGTCTACTTTGGGAAGTTATTAAAGATAATGTTATTCATGAATTGTATCCTAGAAATAAACATGGCAAATATCCTTACTACAGAAGAGCAAGAGACGGAACATGGGGACAATCTAGCCAATTGACTCAGGATGAAGACTTGTGGTTGACTCTAGAGTTTAATCCAGAATTGATGGGGTATGTGTAGGAAGGAGTAAAGAATGGAAACTGGAATTTGTATACGATGTGCTTTCAACTGGGTTACACCTTGCATCAATGAAA
""")

    @patch('bin.get_insertion_sequences.process_sam_file')
    def test_main_diff_min_max_length(self, mock_process_sam_file):
        args = get_arguments().parse_args(
            ['--contig_fasta', self.TEST_CONTIG_FASTA, '--contig_info', self.TEST_CONTIG_FILE, '--sam_file', 'sam_file', '--clstr_file', self.TEST_CLSTR, '--min_is_len', '600', '--max_is_len', '3000', '--output_prefix', self.TEST_OUTPUT_PREFIX_MAIN2])

        mock_process_sam_file.return_value = mapping_info = {
'NODE_1_length_49561_cov_5.4235': {'0': [[10, 30], [550,570], [1,20], [580,600]]},
'NODE_5_length_32495_cov_4.56202': {'8601': [[51, 150]]}}

        main(args)

        txt = open(f'{self.TEST_OUTPUT_PREFIX_MAIN2}_insertion_sequences_info.txt', "r")
        content = "".join(txt.readlines())
        os.remove(f'{self.TEST_OUTPUT_PREFIX_MAIN2}_insertion_sequences_info.txt')
        self.assertEqual(content,"""IS_name	sample_id	contig	is_start	is_end	description
IS_length_600_IPR036515_49_549_IPR036515_52_534_IPR002686_67_522_IPR002686_82_513	tests/data/output/test3	NODE_1_length_49561_cov_5.4235	1	600	Transposase IS200-like superfamily;Transposase IS200-like superfamily;Transposase IS200-like;Transposase IS200-like
""")

        fasta = open(f'{self.TEST_OUTPUT_PREFIX_MAIN2}_insertion_sequences.fasta', "r")
        content = "".join(fasta.readlines())
        os.remove(f'{self.TEST_OUTPUT_PREFIX_MAIN2}_insertion_sequences.fasta')
        self.assertEqual(content,""">IS_length_600_IPR036515_49_549_IPR036515_52_534_IPR002686_67_522_IPR002686_82_513
CAGCAGCCAGTACGGTCTAGGAAGCACCTTTACCCTTGTTCTCAATCTCTCTGGTAGTGAAAATAAAGCTTGTAAGTGTCTAGTGTTTTGAGGTAAATAGAAAATATTGAGCTAACGGAAAATTCCGTTAGTTCTTTTTTGTATTTTAGGAGAAAAAATGTCAGATAATATTGAAATTGTTATGGATATTGATTCGTTAGAAATTTTCGATGTAATAGATCTTGATAAGAATGAAAATATGGTTGATTGTATCCCTGAAGGAAATTCTGAAGAAAAGGTATCTCCATACGATGACCTCAATAGGGTTTGGCTTAGAAACGGGTATTATTATAAGCCGTCTAGTCTACTTTGGGAAGTTATTAAAGATAATGTTATTCATGAATTGTATCCTAGAAATAAACATGGCAAATATCCTTACTACAGAAGAGCAAGAGACGGAACATGGGGACAATCTAGCCAATTGACTCAGGATGAAGACTTGTGGTTGACTCTAGAGTTTAATCCAGAATTGATGGGGTATGTGTAGGAAGGAGTAAAGAATGGAAACTGGAATTTGTATACGATGTGCTTTCAACTGGGTTACACCTTGCATCAATGAAA
""")
