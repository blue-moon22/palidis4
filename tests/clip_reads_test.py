import argparse, os
import unittest
from unittest.mock import patch, call, ANY

from bin.clip_reads import *

class TestClipReads(unittest.TestCase):
    TEST_FASTA = "tests/data/input/test_IR_1.fasta"
    TEST_OUTPUT_PREFIX = "tests/data/output/test"

    def test_clip_reads(self):
        actual = clip_reads(self.TEST_FASTA, self.TEST_OUTPUT_PREFIX)

        fasta_name = self.TEST_OUTPUT_PREFIX + '_irs.fasta'
        fasta = open(fasta_name, "r")
        actual = "".join(fasta.readlines())
        os.remove(fasta_name)
        self.assertEqual(actual, """>test_nstart_test_nend_SRR062032.2555.1_f1_LCoord_29_RCoord_56\nACTATTTATTTTATACCCGACATAAAAT\n""")

    def test_arguments(self):
        actual = get_arguments().parse_args(
            ['--read_fasta', 'read_fasta', '--output_prefix', 'output_prefix'])
        self.assertEqual(actual, argparse.Namespace(read_fasta='read_fasta', output_prefix='output_prefix'))

    @patch('bin.clip_reads.clip_reads')
    def test_main(self, mock_clip_reads):
        args = get_arguments().parse_args(
            ['--read_fasta', 'read_fasta', '--output_prefix', 'output_prefix'])

        main(args)

        self.assertEqual(mock_clip_reads.call_args_list, [call('read_fasta', 'output_prefix')])
