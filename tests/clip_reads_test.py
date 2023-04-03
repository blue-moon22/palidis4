import argparse, os
import unittest
from unittest.mock import patch, call, ANY

from bin.clip_reads import *

class TestClipReads(unittest.TestCase):
    TEST_FASTA = "tests/data/input/test_IR_1.fasta"
    TEST_TAB = "tests/data/input/test_IR.tab"
    TEST_OUTPUT_PREFIX = "tests/data/output/test"

    def test_get_positions(self):
        actual = get_positions(self.TEST_TAB)

        self.assertEqual(actual["Seq142_nstart_SRS011144_nend_SRR062032.4389.1_f1"], [(34,59),(33,59)])

    def test_write_fasta(self):

        positions = get_positions(self.TEST_TAB)
        seq_list = unique_seq(self.TEST_FASTA, positions)
        actual = write_fasta(seq_list, self.TEST_OUTPUT_PREFIX)

        fasta_name = self.TEST_OUTPUT_PREFIX + '_irs.fasta'
        fasta = open(fasta_name, "r")
        actual = "".join(fasta.readlines())
        os.remove(fasta_name)
        self.assertEqual(actual, """>1
TTATTTTATACCCGACATAAAATATC
>2
TTATTTTATACCCGACATAAAATATC
>3
TTTATTTTATACCCGACATAAAATATC
>4
TTTATTTTATACCCGACATAAAATATC
""")

    def test_arguments(self):
        actual = get_arguments().parse_args(
            ['--read_fasta', 'read_fasta', '--tab_file', 'tab_file', '--output_prefix', 'output_prefix'])
        self.assertEqual(actual, argparse.Namespace(read_fasta='read_fasta', tab_file='tab_file', output_prefix='output_prefix'))

    @patch('bin.clip_reads.write_fasta')
    @patch('bin.clip_reads.unique_seq')
    @patch('bin.clip_reads.get_positions')
    def test_main(self, mock_get_positions, mock_unique_seq, mock_write_fasta):
        args = get_arguments().parse_args(
            ['--read_fasta', 'read_fasta', '--tab_file', 'tab_file', '--output_prefix', 'output_prefix'])
        mock_get_positions.return_value = {'a dictionary'}
        mock_unique_seq.return_value = ['a list']

        main(args)

        self.assertEqual(mock_get_positions.call_args_list, [call('tab_file')])
        self.assertEqual(mock_unique_seq.call_args_list, [call('read_fasta', {'a dictionary'})])
        self.assertEqual(mock_write_fasta.call_args_list, [call(['a list'], 'output_prefix')])
