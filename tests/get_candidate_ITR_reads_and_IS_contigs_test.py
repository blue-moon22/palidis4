import argparse, os
import unittest
from unittest.mock import patch, call, ANY

from bin.get_candidate_ITR_reads_and_IS_contigs import allocate_seqs, get_array_size, get_ir_offset, get_contigs_with_itrs, write_fasta_files, write_tab_file, write_contig_file, get_arguments, main

class TestCandidateITRReadsAndISContigs(unittest.TestCase):
    TEST_SAM_FILE1 = 'tests/data/input/test1.sam.mapped.sorted'
    TEST_SAM_FILE2 = 'tests/data/input/test2.sam.mapped.sorted'
    TEST_FASTA1 = 'tests/data/input/test_IR_1.fasta'
    TEST_FASTA2 = 'tests/data/input/test_IR_2.fasta'
    TEST_CONTIG = 'tests/data/input/test.fa'
    TEST_OUTPUT_PREFIX = 'tests/data/output/test'

    size = get_array_size(TEST_FASTA1, TEST_FASTA2)
    offset_arrays = get_ir_offset(TEST_FASTA1, TEST_FASTA2, size)

    def test_get_ir_offset(self):
        actual = get_ir_offset(self.TEST_FASTA1, self.TEST_FASTA2, self.size)
        self.assertEqual(actual[0][48970], 32)

    def test_contigs_with_itrs_in_1(self):
        """
        Test contigs with a pair of ITRs
        """
        alloc_arrays = allocate_seqs(self.size)
        actual = get_contigs_with_itrs(self.TEST_SAM_FILE1, alloc_arrays[0], self.offset_arrays[0], 1, 200.01, 100.01)

        self.assertEqual(actual[0].get('b3_06_1'), None)

    def test_contigs_with_itrs_in_2(self):
        """
        Test contigs with more than two positions within range
        """
        alloc_arrays = allocate_seqs(self.size)
        actual = get_contigs_with_itrs(self.TEST_SAM_FILE2, alloc_arrays[1], self.offset_arrays[1], 2, 200.01, 100.01)

        self.assertEqual(actual[0].get('b3_06_1'), None)

    def test_write_fasta_files(self):
        """
        Test fasta file outputs and returned dictionary
        """
        alloc_arrays = allocate_seqs(self.size)
        output1 = get_contigs_with_itrs(self.TEST_SAM_FILE1, alloc_arrays[0], self.offset_arrays[0], 1, 200.01, 100.01)
        output2 = get_contigs_with_itrs(self.TEST_SAM_FILE2, alloc_arrays[1], self.offset_arrays[1], 2, 200.01, 100.01)
        write_fasta_files(self.TEST_FASTA1, self.TEST_FASTA2, (output1[1], output2[1]), self.TEST_OUTPUT_PREFIX)

        f1_name = self.TEST_OUTPUT_PREFIX + '_reads_with_candidate_itrs_1.fasta'
        f1 = open(f1_name, "r")
        actual = "".join(f1.readlines())
        os.remove(f1_name)
        self.assertEqual(actual, """""")

        f2_name = self.TEST_OUTPUT_PREFIX + '_reads_with_candidate_itrs_2.fasta'
        f2 = open(f2_name, "r")
        actual = "".join(f2.readlines())
        os.remove(f2_name)
        self.assertEqual(actual, """""")

    def test_write_tab_file(self):
        """
        Test tab file output
        """
        alloc_arrays = allocate_seqs(self.size)
        output1 = get_contigs_with_itrs(self.TEST_SAM_FILE1, alloc_arrays[0], self.offset_arrays[0], 1, 200.01, 100.01)
        output2 = get_contigs_with_itrs(self.TEST_SAM_FILE2, alloc_arrays[1], self.offset_arrays[1], 2, 200.01, 100.01)
        write_tab_file(output1[0], output2[0], self.TEST_OUTPUT_PREFIX)

        tab_name = self.TEST_OUTPUT_PREFIX + '_contigs_reads_ir_position_info.tab'
        tab = open(tab_name, "r")
        actual = "".join(tab.readlines())
        os.remove(tab_name)
        self.assertEqual(actual, """sample_id\tcontig\tread\tposition\n""")

    def test_write_contig_file(self):
        """
        Test contig file output
        """
        itr_info1 = {'b3_06_1': {'Seq8080_nstart_b3_06_nend_HS38_13470:2:1101:7938:91791#84/1_f1_LCoord_91_RCoord_116': 1776263, 'Seq48971_nstart_b3_06_nend_HS38_13470:2:1102:19036:87556#84/2_f2': 1775404}}
        itr_info2 = {'b3_06_1': {'Seq8080_nstart_b3_06_nend_HS38_13470:2:1101:7938:91791#84/1_f1_LCoord_91_RCoord_116': 1776263, 'Seq48971_nstart_b3_06_nend_HS38_13470:2:1102:19036:87556#84/2_f2': 1775404}}
        write_contig_file(self.TEST_CONTIG, itr_info1, itr_info2, self.TEST_OUTPUT_PREFIX)
        contig_file_name = self.TEST_OUTPUT_PREFIX + '_contigs_with_candidate_itrs.fa'
        file = open(contig_file_name, "r")
        actual = "".join(file.readlines())
        os.remove(contig_file_name)
        self.assertEqual(actual, """>b3_06_1\nTTGGTGGATTTGTTGGTGGATTTGTTGGGGGATTTGTTGGTGGATTTGTTGGTGGATTTG\nTTGGTGGATTTGTTGGTGCACTATTTACTACCACTTTCACTTTAACATTCGCAAGTTTAC\nAATCAGTACCGTACTTTGTGCAAAGCTCATATTTAAATCGATAACTTCCTTGCGGATATT\n""")

    def test_arguments(self):
        actual = get_arguments().parse_args(
            ['--contig_fasta', 'contig_path', '--sam_file1', 'sam_path1', '--sam_file2', 'sam_path2', '--fasta1', 'fasta1_path',
            '--fasta2', 'fasta2_path', '--insert_size', '200.0', '--read_length', '100.0', '--output_prefix', 'output_prefix'])
        self.assertEqual(actual,
                         argparse.Namespace(contig_fasta='contig_path', sam_file1='sam_path1', sam_file2='sam_path2', fasta_file1='fasta1_path', fasta_file2='fasta2_path', insert_size=200.0, read_length=100.0, output_prefix='output_prefix'))

    def test_arguments_short_options(self):
        actual = get_arguments().parse_args(
            ['-c', 'contig_path', '-s1', 'sam_path1', '-s2', 'sam_path2', '-f1', 'fasta1_path',
            '-f2', 'fasta2_path', '-i', '200.0', '-l', '100.0', '-o', 'output_prefix'])
        self.assertEqual(actual,
                         argparse.Namespace(contig_fasta='contig_path', sam_file1='sam_path1', sam_file2='sam_path2', fasta_file1='fasta1_path', fasta_file2='fasta2_path', insert_size=200.0, read_length=100.0, output_prefix='output_prefix'))

@patch('bin.get_candidate_ITR_reads_and_IS_contigs.get_array_size')
@patch('bin.get_candidate_ITR_reads_and_IS_contigs.allocate_seqs')
@patch('bin.get_candidate_ITR_reads_and_IS_contigs.get_ir_offset')
@patch('bin.get_candidate_ITR_reads_and_IS_contigs.get_contigs_with_itrs')
@patch('bin.get_candidate_ITR_reads_and_IS_contigs.get_contigs_with_itrs')
@patch('bin.get_candidate_ITR_reads_and_IS_contigs.write_contig_file')
@patch('bin.get_candidate_ITR_reads_and_IS_contigs.write_fasta_files')
@patch('bin.get_candidate_ITR_reads_and_IS_contigs.write_tab_file')
def test_main(self, mock_write_tab_file, mock_write_fasta_files, mock_write_contig_file, mock_get_contigs_with_itrs, mock_get_ir_offset, mock_allocate_seqs, mock_get_array_size):
    args = get_arguments().parse_args(
        ['--contig_fasta', 'contig_path', '--sam_file1', 'sam_path1', '--sam_file2', 'sam_path2', '--fasta1', 'fasta1_path',
        '--fasta2', 'fasta2_path', '--insert_size', '200.0', '--read_length', '100.0', '--output_prefix', 'output_prefix'])
    mock_get_array_size.return_value = 1000
    mock_allocate_seqs.return_value = ['alloc1', 'alloc2']
    mock_get_ir_offset.return_value = ['offset1', 'offset2']
    mock_get_contigs_with_itrs.side_effect = [[{'contigs': {'reads': 'positions'}}, ['alloc']], [{'contigs': {'reads': 'positions'}}, ['alloc']]]

    main(args)

    mock_get_array_size.call_args_list = call('fasta1_path', 'fasta2_path')
    mock_allocate_seqs.call_args_list = 1000
    mock_get_ir_offset.call_args_list = call('fasta1_path', 'fasta2_path', 1000)
    mock_get_contigs_with_itrs.call_count = 2
    mock_get_contigs_with_itrs.call_args_list = [call('sam_path1', 'alloc1', 'offset1', 1, 200.0, 100.0), call('sam_path2', 'alloc2', 'offset2', 2, 200.0, 100.0)]
    mock_write_contig_file.call_args_list = call('contig_path', {'contigs': {'reads': 'positions'}}, {'contigs': {'reads': 'positions'}}, 'output_prefix')
    mock_write_fasta_files.call_args_list = call('fasta1_path', 'fasta2_path', (['alloc'], ['alloc']), 'output_prefix')
    mock_write_tab_file.call_args_list = call({'contigs': {'reads': 'positions'}}, {'contigs': {'reads': 'positions'}}, 'output_prefix')
