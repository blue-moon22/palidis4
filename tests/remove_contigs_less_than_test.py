import argparse, os
import unittest
from unittest.mock import patch, call, ANY

from bin.remove_contigs_smaller_than import *

class TestRemoveContigsLessThanTest(unittest.TestCase):
    TEST_CONTIG_FILE = 'tests/data/input/test_scaffolds.fasta'
    TEST_OUTPUT_FILE = 'tests/data/output/test_scaffolds_filtered.fasta'

    def test_select_contigs(self):
        """
        Remove contigs
        """
        select_contigs(self.TEST_CONTIG_FILE, 500, self.TEST_OUTPUT_FILE)

        file = open(self.TEST_OUTPUT_FILE, "r")
        actual = "".join(file.readlines())
        os.remove(self.TEST_OUTPUT_FILE)

        self.assertEqual(actual.count('>'), 7281)


    def test_arguments(self):
        actual = get_arguments().parse_args(
            ['--min_length', '500', '--fasta_file', self.TEST_CONTIG_FILE, '--output_file', self.TEST_OUTPUT_FILE])
        self.assertEqual(actual,
                         argparse.Namespace(fasta_file=self.TEST_CONTIG_FILE, min_length=500, output_file=self.TEST_OUTPUT_FILE))


    @patch('bin.remove_contigs_smaller_than')
    def test_main(self, mock_select_contigs):
        args = get_arguments().parse_args(['--min_length', '500', '--fasta_file', self.TEST_CONTIG_FILE, '--output_file', self.TEST_OUTPUT_FILE])

        main(args)

        mock_select_contigs.call_args_list = call(self.TEST_CONTIG_FILE, 500, self.TEST_OUTPUT_FILE)
