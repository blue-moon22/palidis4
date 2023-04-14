#!/usr/bin/env python3

"""
author: Victoria Carr

Get candidate insertion sequences.
"""

import argparse, sys
import re
from collections import defaultdict
import logging

class contigTransposaseInfo:

    def __init__(self):
        self.contig_info = {}

    def load_contig_info(self, contig_info_file):
        with open(contig_info_file, "r") as file:
            next(file)
            for line in file:
                protein = line.split("\t")[0]
                accession = line.split("\t")[1]
                annotation = line.split("\t")[2]
                start = int(line.split("\t")[3])
                end = int(line.split("\t")[4].replace("\n", ""))
                contig = protein[:protein.rindex('_')]

                if contig not in self.contig_info:
                    self.contig_info[contig] = []
                tmp = self.contig_info[contig]
                tmp.append((accession, annotation, start, end))
                self.contig_info[contig] = tmp

    def get_contig_info(self):
        return(self.contig_info)

    def check_irs_flanking(self, contig, ir_positions, MIN_IS_LEN, MAX_IS_LEN):
        is_info = {}
        for index1, ir_positions1 in enumerate(ir_positions[0]):
            for index2, ir_positions2 in enumerate(ir_positions[1]):
                if ir_positions2[0] > ir_positions1[0]:
                    length = (ir_positions2[1] + 1) - ir_positions1[0]
                    itr1_start = ir_positions1[0]
                    itr1_end = ir_positions1[1]
                    itr2_start = ir_positions2[0]
                    itr2_end = ir_positions2[1]
                else:
                    length = (ir_positions1[1] + 1) - ir_positions2[0]
                    itr1_start = ir_positions2[0]
                    itr1_end = ir_positions2[1]
                    itr2_start = ir_positions1[0]
                    itr2_end = ir_positions1[1]

                if length >= MIN_IS_LEN and length <= MAX_IS_LEN:
                    for transposase in self.contig_info[contig]:
                        start = transposase[2]
                        end = transposase[3]
                        if start > itr1_end and end < itr2_start:
                            if (itr1_start, itr2_end) not in is_info:
                                is_info[(itr1_start, itr2_end)] = []
                            tmp = is_info[(itr1_start, itr2_end)]
                            tmp.append(transposase)
                            tmp = list(set(tmp))
                            tmp.sort()
                            is_info[(itr1_start, itr2_end)] = tmp

        return is_info

def reverse_complement(dna_sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    reverse_complement_seq = ''
    for base in dna_sequence[::-1]:
        reverse_complement_seq += complement_dict[base]
    return reverse_complement_seq

def get_candidate_itrs(sam_file):

    with open(sam_file, "r") as sam:

        itr_candidates = {}
        for line in sam:
            sequence = line.split("\t")[9]
            contig = line.split("\t")[2]
            rev_comp = reverse_complement(sequence)

            if (contig,sequence) not in itr_candidates:
                itr_candidates[(contig,sequence)] = False

            if (contig,rev_comp) in itr_candidates:
                itr_candidates[(contig,rev_comp)] = True

        irs = []
        for key, item in itr_candidates.items():
            if item:
                rev_comp = reverse_complement(key[1])
                irs.append(key)
                irs.append((key[0], rev_comp))

        irs = list(set(irs))
        irs.sort()

    return irs


def process_sam_file(sam_file, irs):
    """
    Function to read information from the sam file and assign positions of IRs
    """

    with open(sam_file, "r") as sam:

        mapping_info = {}
        for line in sam:
            contig = line.split("\t")[2]
            sequence = line.split("\t")[9]

            if (contig,sequence) in irs:

                contig = line.split("\t")[2]
                pos = int(line.split("\t")[3])
                len_fragment = len(line.split("\t")[9])
                ir_positions = (pos, (pos + len_fragment)-1)

                ir_pairs = [sequence, reverse_complement(sequence)]
                ir_pairs.sort()
                ir = ir_pairs[0]
                if contig not in mapping_info:
                    mapping_info[contig] = {}
                    mapping_info[contig][ir] = [[],[]]
                else:
                    if ir not in mapping_info[contig]:
                        mapping_info[contig][ir] = [[],[]]

                tmp = mapping_info[contig][ir]
                if sequence == ir:
                    if ir_positions not in tmp[0]:
                        tmp[0].append(ir_positions)
                else:
                    if ir_positions not in tmp[1]:
                        tmp[1].append(ir_positions)
                mapping_info[contig][ir] = tmp

    return mapping_info


def check_itrs(contig, ir_positions, MIN_IS_LEN, MAX_IS_LEN):
    is_info = []
    for index1, ir_positions1 in enumerate(ir_positions[0]):
        for index2, ir_positions2 in enumerate(ir_positions[1]):
            if ir_positions2[0] > ir_positions1[0]:
                length = (ir_positions2[1] + 1) - ir_positions1[0]
                itr1_start = ir_positions1[0]
                itr1_end = ir_positions1[1]
                itr2_start = ir_positions2[0]
                itr2_end = ir_positions2[1]
            else:
                length = (ir_positions1[1] + 1) - ir_positions2[0]
                itr1_start = ir_positions2[0]
                itr1_end = ir_positions2[1]
                itr2_start = ir_positions1[0]
                itr2_end = ir_positions1[1]

            if length >= MIN_IS_LEN and length <= MAX_IS_LEN:
                is_info.append((itr1_start, itr2_end))

    return is_info


def write_insertion_sequences_info(mapping_info, sample_id, MIN_IS_LEN, MAX_IS_LEN):

    contig_is_info = {}
    with open(f'{sample_id}_candidate_insertion_sequences_info.txt', "w") as out:
        out.write("sample_id\tcontig\titr1_start\titr2_end\n")
        for contig, ir in mapping_info.items():
            count = 1
            for seq, ir_positions in ir.items():
                itrs = check_itrs(contig, ir_positions, MIN_IS_LEN, MAX_IS_LEN)
                for itr in itrs:
                    contig_name = f'{contig}_{str(count)}'
                    out.write(f'{sample_id}\t{contig_name}\t{itr[0]}\t{itr[1]}\n')
                    if contig not in contig_is_info:
                        contig_is_info[contig] = []
                    tmp = contig_is_info[contig]
                    itr = list(itr)
                    itr.append(contig_name)
                    tmp.append(itr)
                    contig_is_info[contig] = tmp
                    count += 1

    return contig_is_info


def write_insertion_sequences_fasta(contig_is_info, contig_fasta, output_prefix):

    contig = None
    with open(f'{output_prefix}_candidate_insertion_sequences.fasta', "w") as out:
        with open(contig_fasta, "r") as fasta:
            for line in fasta:
                if line[0] == ">":
                    contig = line.replace(">", "").replace("\n", "")
                else:
                    if contig in contig_is_info:
                        for itrs in contig_is_info[contig]:
                            print(itrs)
                            out.write(f'>{itrs[2]}\n{line[itrs[0]-1:(itrs[1])]}\n')


def get_arguments():

    parser = argparse.ArgumentParser(description='Get reads that contain candidate ITRs and contigs that associate with these.')
    parser.add_argument('--contig_fasta', '-f', dest='contig_fasta', required=True,
                        help='Input contig FASTA file.', type = str)
    parser.add_argument('--sam_file', '-s', dest='sam_file', required=True,
                        help='Input sam file.', type = str)
    parser.add_argument('--min_is_len', '-min', dest='min_is_len', required=False,
                        help='Minimum length of insertion sequence', type = int, default = 500)
    parser.add_argument('--max_is_len', '-max', dest='max_is_len', required=False,
                        help='Maximum length of insertion sequence', type = int, default = 3000)
    parser.add_argument('--output_prefix', '-o', dest='output_prefix', required=True,
                        help='Prefix of output files (should be sample ID).', type = str)
    return parser


def main(args):

    logging.basicConfig(
            format='%(asctime)s %(levelname)-8s %(message)s',
            level=logging.INFO,
            datefmt='%Y-%m-%d %H:%M:%S')

    # Get IR mapping info
    logging.info('Get candidate ITRs')
    irs = get_candidate_itrs(args.sam_file)
    logging.info('Process SAM file')
    mapping_info = process_sam_file(args.sam_file, irs)

    # Write insertion sequence info
    logging.info('Write candidate insertion sequence info')
    contig_is_info = write_insertion_sequences_info(mapping_info, args.output_prefix, args.min_is_len, args.max_is_len)

    # Write insertion sequence fasta
    print(contig_is_info)
    logging.info('Write FASTA file')
    write_insertion_sequences_fasta(contig_is_info, args.contig_fasta, args.output_prefix)


if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))
