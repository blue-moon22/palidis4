import argparse, os
import unittest
from unittest.mock import patch, call, ANY

from bin.get_read_and_fasta_info import *

class GetReadAndFastaInfo(unittest.TestCase):
    TEST_ANNOT = 'tests/data/input/test_insertion_sequence_annotations.tab'
    TEST_IR_TAB = 'tests/data/input/test_reads_itr_clusters.txt'
    TEST_FASTA = 'tests/data/input/test2_irs.fasta'
    TEST_OUTPUT_PREFIX = 'tests/data/output/test'

    headers = ['Seq377725_nstart_SRS012281_nend_SRR059448.31040150.1_f1_LCoord_29_RCoord_54',
        'Seq442389_nstart_SRS012281_nend_SRR059449.1386939.1_f1_LCoord_48_RCoord_75',
        'Seq499624_nstart_SRS012281_nend_SRR059449.5536288.1_f1_LCoord_35_RCoord_62',
        'Seq77063_nstart_SRS012281_nend_SRR059448.6504893.1_f1_LCoord_31_RCoord_74',
        'Seq362409_nstart_SRS012281_nend_SRR059448.29723236.1_f1_LCoord_39_RCoord_73',
        'Seq337950_nstart_SRS012281_nend_SRR059448.27648769.1_f1_LCoord_49_RCoord_74',
        'Seq860915_nstart_SRS012281_nend_SRR059449.31901904.2_f2_LCoord_41_RCoord_73',
        'Seq195412_nstart_SRS012281_nend_SRR059448.16075499.2_f2_LCoord_51_RCoord_76',
        'Seq484253_nstart_SRS012281_nend_SRR059449.4428088.2_f2_LCoord_51_RCoord_76',
        'Seq722376_nstart_SRS012281_nend_SRR059449.21595447.2_f2_LCoord_35_RCoord_62',
        'Seq652277_nstart_SRS012281_nend_SRR059449.16517467.2_f2_LCoord_31_RCoord_75',
        'Seq18815_nstart_SRS012281_nend_SRR059448.1573381.2_f2_LCoord_49_RCoord_74',
        'Seq651223_nstart_SRS012281_nend_SRR059449.16445107.2_f2_LCoord_50_RCoord_75',
        'Seq803343_nstart_SRS012281_nend_SRR059449.27647569.1_f1_LCoord_28_RCoord_67',
        'Seq229631_nstart_SRS012281_nend_SRR059448.18806974.2_f2_LCoord_29_RCoord_68',
        'Seq151365_nstart_SRS012281_nend_SRR059448.12500609.2_f2_LCoord_40_RCoord_74',
        'Seq747627_nstart_SRS012281_nend_SRR059449.23488848.1_f1_LCoord_45_RCoord_76',
        'Seq243173_nstart_SRS012281_nend_SRR059448.19882929.2_f2_LCoord_32_RCoord_59',
        'Seq734135_nstart_SRS012281_nend_SRR059449.22490871.2_f2_LCoord_31_RCoord_58']

    def test_get_itr_clusters(self):
        actual = get_itr_clusters(self.TEST_ANNOT)

        self.assertEqual(actual, ['NODE_7_length_63398_cov_8.5468549','NODE_30_length_27264_cov_6.38506210','NODE_49_length_14427_cov_6.96277250'])


    def test_get_read_headers_of_clusters(self):
        actual = get_read_headers_of_clusters(['NODE_7_length_63398_cov_8.5468549','NODE_30_length_27264_cov_6.38506210','NODE_49_length_14427_cov_6.96277250'], self.TEST_IR_TAB, self.TEST_OUTPUT_PREFIX)

        self.assertEqual(actual, self.headers)

        file_name = self.TEST_OUTPUT_PREFIX + '_itr_read_positions_clusters.txt'
        file = open(file_name, "r")
        actual = "".join(file.readlines())
        os.remove(file_name)
        self.maxDiff = None
        self.assertEqual(actual, """sample_id\tcontig\tread\titr_start\titr_end\titr_cluster\nSRS012281\tNODE_7_length_63398_cov_8.54685\tSeq377725_nstart_SRS012281_nend_SRR059448.31040150.1_f1_LCoord_29_RCoord_54\t27516\t27541\t49\nSRS012281\tNODE_7_length_63398_cov_8.54685\tSeq442389_nstart_SRS012281_nend_SRR059449.1386939.1_f1_LCoord_48_RCoord_75\t31195\t31222\t49\nSRS012281\tNODE_7_length_63398_cov_8.54685\tSeq499624_nstart_SRS012281_nend_SRR059449.5536288.1_f1_LCoord_35_RCoord_62\t31195\t31222\t49\nSRS012281\tNODE_7_length_63398_cov_8.54685\tSeq77063_nstart_SRS012281_nend_SRR059448.6504893.1_f1_LCoord_31_RCoord_74\t56515\t56558\t49\nSRS012281\tNODE_7_length_63398_cov_8.54685\tSeq362409_nstart_SRS012281_nend_SRR059448.29723236.1_f1_LCoord_39_RCoord_73\t56527\t56561\t49\nSRS012281\tNODE_7_length_63398_cov_8.54685\tSeq337950_nstart_SRS012281_nend_SRR059448.27648769.1_f1_LCoord_49_RCoord_74\t56533\t56558\t49\nSRS012281\tNODE_7_length_63398_cov_8.54685\tSeq860915_nstart_SRS012281_nend_SRR059449.31901904.2_f2_LCoord_41_RCoord_73\t18277\t18309\t49\nSRS012281\tNODE_7_length_63398_cov_8.54685\tSeq195412_nstart_SRS012281_nend_SRR059448.16075499.2_f2_LCoord_51_RCoord_76\t29993\t30018\t49\nSRS012281\tNODE_7_length_63398_cov_8.54685\tSeq484253_nstart_SRS012281_nend_SRR059449.4428088.2_f2_LCoord_51_RCoord_76\t29994\t30019\t49\nSRS012281\tNODE_7_length_63398_cov_8.54685\tSeq722376_nstart_SRS012281_nend_SRR059449.21595447.2_f2_LCoord_35_RCoord_62\t31195\t31222\t49\nSRS012281\tNODE_7_length_63398_cov_8.54685\tSeq652277_nstart_SRS012281_nend_SRR059449.16517467.2_f2_LCoord_31_RCoord_75\t50516\t50560\t49\nSRS012281\tNODE_7_length_63398_cov_8.54685\tSeq18815_nstart_SRS012281_nend_SRR059448.1573381.2_f2_LCoord_49_RCoord_74\t56533\t56558\t49\nSRS012281\tNODE_7_length_63398_cov_8.54685\tSeq651223_nstart_SRS012281_nend_SRR059449.16445107.2_f2_LCoord_50_RCoord_75\t56533\t56558\t49\nSRS012281\tNODE_30_length_27264_cov_6.38506\tSeq803343_nstart_SRS012281_nend_SRR059449.27647569.1_f1_LCoord_28_RCoord_67\t14425\t14464\t210\nSRS012281\tNODE_30_length_27264_cov_6.38506\tSeq229631_nstart_SRS012281_nend_SRR059448.18806974.2_f2_LCoord_29_RCoord_68\t14425\t14464\t210\nSRS012281\tNODE_30_length_27264_cov_6.38506\tSeq151365_nstart_SRS012281_nend_SRR059448.12500609.2_f2_LCoord_40_RCoord_74\t17079\t17113\t210\nSRS012281\tNODE_49_length_14427_cov_6.96277\tSeq747627_nstart_SRS012281_nend_SRR059449.23488848.1_f1_LCoord_45_RCoord_76\t2661\t2692\t250\nSRS012281\tNODE_49_length_14427_cov_6.96277\tSeq243173_nstart_SRS012281_nend_SRR059448.19882929.2_f2_LCoord_32_RCoord_59\t18\t45\t250\nSRS012281\tNODE_49_length_14427_cov_6.96277\tSeq734135_nstart_SRS012281_nend_SRR059449.22490871.2_f2_LCoord_31_RCoord_58\t18\t45\t250\n""")

    def test_write_itr_reads(self):

        write_itr_reads(self.headers, self.TEST_FASTA, self.TEST_OUTPUT_PREFIX)

        fasta_name = self.TEST_OUTPUT_PREFIX + '_ITRs.fasta'
        fasta = open(fasta_name, "r")
        actual = "".join(fasta.readlines())
        os.remove(fasta_name)
        self.maxDiff = None
        self.assertEqual(actual, """>Seq77063_nstart_SRS012281_nend_SRR059448.6504893.1_f1_LCoord_31_RCoord_74\nTACTCAATGAAAATCAAAGAGCAAACTAGGAAACTAGCCACAGG\n>Seq337950_nstart_SRS012281_nend_SRR059448.27648769.1_f1_LCoord_49_RCoord_74\nTACTCAATGAAAATCAAAGAGCAAAC\n>Seq362409_nstart_SRS012281_nend_SRR059448.29723236.1_f1_LCoord_39_RCoord_73\nTTATACTCAATGAAAATCAAAGAGCAAACTAGGAA\n>Seq377725_nstart_SRS012281_nend_SRR059448.31040150.1_f1_LCoord_29_RCoord_54\nTTTGCTCTTTGATTTTCATTGAGTAT\n>Seq442389_nstart_SRS012281_nend_SRR059449.1386939.1_f1_LCoord_48_RCoord_75\nTTATACTCAATGAAAATCAAAGAGCAAA\n>Seq499624_nstart_SRS012281_nend_SRR059449.5536288.1_f1_LCoord_35_RCoord_62\nTTATACTCAATGAAAATCAAAGAGCAAA\n>Seq747627_nstart_SRS012281_nend_SRR059449.23488848.1_f1_LCoord_45_RCoord_76\nATGAAAATCAAAGAGCAAACTAGGAAGCTAGC\n>Seq803343_nstart_SRS012281_nend_SRR059449.27647569.1_f1_LCoord_28_RCoord_67\nCAACCTCAAAACAGTGTTTTGAGCTGACTTCGTCAGTCTT\n>Seq18815_nstart_SRS012281_nend_SRR059448.1573381.2_f2_LCoord_49_RCoord_74\nTACTCAATGAAAATCAAAGAGCAAAC\n>Seq151365_nstart_SRS012281_nend_SRR059448.12500609.2_f2_LCoord_40_RCoord_74\nTCAAAACAGTGTTTTGAGCTGACTTCGTCAGTCTT\n>Seq195412_nstart_SRS012281_nend_SRR059448.16075499.2_f2_LCoord_51_RCoord_76\nGTTTGCTCTTTGATTTTCATTGAGTA\n>Seq229631_nstart_SRS012281_nend_SRR059448.18806974.2_f2_LCoord_29_RCoord_68\nCAACCTCAAAACAGTGTTTTGAGCTGACTTCGTCAGTCTT\n>Seq243173_nstart_SRS012281_nend_SRR059448.19882929.2_f2_LCoord_32_RCoord_59\nTAGCTTCCTAGTTTGCTCTTTGATTTTC\n>Seq484253_nstart_SRS012281_nend_SRR059449.4428088.2_f2_LCoord_51_RCoord_76\nATACTCAATGAAAATCAAAGAGCAAA\n>Seq651223_nstart_SRS012281_nend_SRR059449.16445107.2_f2_LCoord_50_RCoord_75\nTACTCAATGAAAATCAAAGAGCAAAC\n>Seq652277_nstart_SRS012281_nend_SRR059449.16517467.2_f2_LCoord_31_RCoord_75\nATACTCAATGAAAATCAAAGAGCAAACTAGGAAACTAGCCACAGG\n>Seq722376_nstart_SRS012281_nend_SRR059449.21595447.2_f2_LCoord_35_RCoord_62\nTTATACTCAATGAAAATCAAAGAGCAAA\n>Seq734135_nstart_SRS012281_nend_SRR059449.22490871.2_f2_LCoord_31_RCoord_58\nTAGCTTCCTAGTTTGCTCTTTGATTTTC\n>Seq860915_nstart_SRS012281_nend_SRR059449.31901904.2_f2_LCoord_41_RCoord_73\nATACTCAATGAAAATCAAAGAGCAAACTAGGAA\n""")

    def test_arguments(self):
        actual = get_arguments().parse_args(
            ['--ir_fasta', 'ir_fasta', '--ir_cluster_tab', 'tab_file',
            '--is_annotations_tab', 'is_annot', '--output_prefix', 'output_prefix'])
        self.assertEqual(actual, argparse.Namespace(fasta_file='ir_fasta', tab_file='tab_file', annotations_file='is_annot', output_prefix='output_prefix'))

    def test_main(self):
        args = get_arguments().parse_args(
            ['--ir_fasta', self.TEST_FASTA, '--ir_cluster_tab', self.TEST_IR_TAB,
            '--is_annotations_tab', self.TEST_ANNOT, '--output_prefix', self.TEST_OUTPUT_PREFIX])

        main(args)

        fasta_name = self.TEST_OUTPUT_PREFIX + '_ITRs.fasta'
        os.remove(fasta_name)

        file_name = self.TEST_OUTPUT_PREFIX + '_itr_read_positions_clusters.txt'
        os.remove(file_name)
