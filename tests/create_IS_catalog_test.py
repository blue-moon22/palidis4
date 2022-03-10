import argparse, os
import unittest
from unittest.mock import patch, call, ANY

from bin.create_IS_catalog import *

class TestCreateISCatalog(unittest.TestCase):
    TEST_CLUSTERS = 'tests/data/input/all.fasta.clstr'
    TEST_ANNOT = ['tests/data/input/SRS012281_insertion_sequence_annotations.tab', 'tests/data/input/SRS015040_insertion_sequence_annotations.tab']
    TEST_ITRS = ['tests/data/input/SRS012281_itr_read_positions_clusters.txt', 'tests/data/input/SRS015040_itr_read_positions_clusters.txt']
    TEST_OUTPUT_PREFIX = 'tests/data/output/test'

    def test_get_new_clusters(self):

        actual = get_new_clusters(self.TEST_CLUSTERS)

        self.assertEqual(actual['Seq873926_nstart_SRS017814_nend_SRR060373.2107882.1_f1_LCoord_26_RCoord_75'], '0')
        self.assertEqual(actual['Seq1434744_nstart_SRS017814_nend_SRR060373.3527921.1_f1_LCoord_43_RCoord_74'], '0')
        self.assertEqual(actual['Seq309768_nstart_SRS017127_nend_SRR061945.1520564.1_f1_LCoord_38_RCoord_69'], '1')
        self.assertEqual(actual['Seq3892880_nstart_SRS015040_nend_SRR062525.12815266.1_f1_LCoord_46_RCoord_74'], '53')

    def test_get_itr_reads(self):

        actual = get_itr_reads(self.TEST_ITRS)

        self.maxDiff = None
        self.assertEqual(actual['SRS012281_NODE_7_length_63398_cov_8.54685_49'][0], ('Seq377725_nstart_SRS012281_nend_SRR059448.31040150.1_f1_LCoord_29_RCoord_54', 27516, 27541))


    def test_assign_new_clusters(self):

        assign_new_clusters(get_new_clusters(self.TEST_CLUSTERS), get_itr_reads(self.TEST_ITRS), self.TEST_ANNOT, self.TEST_OUTPUT_PREFIX)

        file_name = self.TEST_OUTPUT_PREFIX + '_insertion_sequence_catalog.txt'
        file = open(file_name, "r")
        actual = "".join(file.readlines())
        os.remove(file_name)
        self.maxDiff = None
        self.assertEqual(actual, """sample_id\tcontig\titr1_start_position\titr1_end_position\titr2_start_position\titr2_end_position\titr_cluster\nSRS012281\tNODE_7_length_63398_cov_8.54685\t29993\t30019\t31195\t31222\t103;56\nSRS012281\tNODE_30_length_27264_cov_6.38506\t14425\t14464\t17079\t17113\t163\nSRS012281\tNODE_49_length_14427_cov_6.96277\t18\t45\t2661\t2692\t56\nSRS015040\tNODE_471_length_6371_cov_62.8583\t1300\t1335\t3955\t3980\t53\nSRS015040\tNODE_471_length_6371_cov_62.8583\t1299\t1324\t3947\t3991\t53\n""")
