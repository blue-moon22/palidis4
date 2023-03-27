import argparse, os
import unittest
from unittest.mock import patch, call, ANY

from bin.convert_fastq_to_fasta import *

class TestConvertFastqToFasta(unittest.TestCase):
    TEST_FASTQ = "tests/data/input/test.fastq"
    TEST_FASTA = "tests/data/output/test.fasta"

    def test_write_fasta(self):

        write_fasta(self.TEST_FASTQ, 1, self.TEST_FASTA)

        fasta_name = self.TEST_FASTA
        fasta = open(fasta_name, "r")
        actual = "".join(fasta.readlines())
        os.remove(fasta_name)
        self.assertEqual(actual, """>Seq1_f1
GTTCCCAACCAGTTGCCATGATGCCTCCATCATTGAAGCGATAGAGTTTACCATCAATTTTATTTAAACCAGTGAACATTTCTCCAGTTTTAGGATCTAGATAATACCATTTACCTGCATCATTGAGCCAACCACGTTTCATTTCACCTGA
>Seq2_f1
ATTGTGGACCTTGTTTTTCTACACGGTTATATTAAGCATCGTGATTTATTCTCTTGTCACAGATTTTTCAAACATCCAAGAATTCATTTACAGTGAATTTTAAATATTTCCTTTAGAAAAACTCCTTGTCTCTAATATTTATTCGAGACA
""")

    def test_arguments(self):
        actual = get_arguments().parse_args(
            ['--fastq', 'fastq_file', '--read', '1', '--fasta', 'fasta_file'])
        self.assertEqual(actual, argparse.Namespace(fastq='fastq_file', read='1', fasta='fasta_file'))
