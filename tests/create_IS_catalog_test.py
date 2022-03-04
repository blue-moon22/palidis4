import argparse, os
import unittest
from unittest.mock import patch, call, ANY

from bin.create_IS_catalog import *

class TestCreateISCatalog(unittest.TestCase):
    TEST_CLIPPED_READS = 'tests/data/input/test_irs.fasta'
    TEST_CLSTR1 = 'tests/data/input/test.fasta.clstr'
    TEST_CLSTR2 = 'tests/data/input/test2.fasta.clstr'
    TEST_INFO_TAB1 = 'tests/data/input/test_contigs_reads_ir_position_info.tab'
    TEST_INFO_TAB2 = 'tests/data/input/test2_contigs_reads_ir_position_info.tab'
    TEST_OUTPUT_PREFIX = 'tests/data/output/test'
    TEST_CONTIG_FASTA = 'tests/data/input/test_assemblies_filtered.fasta'
