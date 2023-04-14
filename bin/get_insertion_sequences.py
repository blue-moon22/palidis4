#!/usr/bin/env python3

"""
author: Victoria Carr

Functions to get insertion sequence information.
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


def write_insertion_sequences_info(contig_transposase_info, contig_info, sample_id):

    contig_is_info = {}
    with open(f'{sample_id}_insertion_sequences_info.txt', "w") as out:
        out.write("IS_name\tsample_id\tcontig\tis_start\tis_end\tdescription\n")
        with open(contig_info, "r") as f:
            for line in f:
                contig = line.split("\t")[1]
                is_start = line.split("\t")[2]
                is_end = line.split("\t")[3].replace("\n", "")

                transp_info = contig_transposase_info.get_contig_info()
                is_name = f'IS_{sample_id}_{contig}_{is_start}_{is_end}'
                description = ''
                if contig in transp_info:
                    for transposase in transp_info[contig]:
                        if description:
                            description = (f'{description};{transposase[0]}:{transposase[1]}:{transposase[2]}:{transposase[3]}')
                        else:
                            description = (f'{transposase[0]}:{transposase[1]}:{transposase[2]}:{transposase[3]}')
                    out.write(f'{is_name}\t{sample_id}\t{contig}\t{is_start}\t{is_end}\t{description}\n')
                    contig_is_info[contig] = is_name

    return contig_is_info


def write_insertion_sequences_fasta(contig_is_info, contig_fasta, output_prefix):

    contig = None
    with open(f'{output_prefix}_insertion_sequences.fasta', "w") as out:
        with open(contig_fasta, "r") as fasta:
            for line in fasta:
                if line[0] == ">":
                    contig = line.replace(">", "").replace("\n", "")
                else:
                    if contig in contig_is_info:
                        out.write(f'>{contig_is_info[contig]}\n{line}')


def get_arguments():

    parser = argparse.ArgumentParser(description='Get reads that contain candidate ITRs and contigs that associate with these.')
    parser.add_argument('--contig_fasta', '-f', dest='contig_fasta', required=True,
                        help='Input contig FASTA file.', type = str)
    parser.add_argument('--contig_info', '-i', dest='contig_info', required=True,
                        help='Input txt file with contig and IS positions.', type = str)
    parser.add_argument('--interproscan_info', '-t', dest='interpro_info', required=True,
                        help='Input tsv file with transposase information.', type = str)
    parser.add_argument('--output_prefix', '-o', dest='output_prefix', required=True,
                        help='Prefix of output files (should be sample ID).', type = str)
    return parser


def main(args):

    logging.basicConfig(
            format='%(asctime)s %(levelname)-8s %(message)s',
            level=logging.INFO,
            datefmt='%Y-%m-%d %H:%M:%S')

    # Get transposase info
    logging.info('Get transpoase info')
    contig_transposase_info = contigTransposaseInfo()
    contig_transposase_info.load_contig_info(args.interpro_info)

    # Write insertion sequence info
    logging.info('Write insertion sequence info')
    contig_is_info = write_insertion_sequences_info(contig_transposase_info, args.contig_info, args.output_prefix)

    # Write insertion sequence fasta
    logging.info('Write FASTA file')
    write_insertion_sequences_fasta(contig_is_info, args.contig_fasta, args.output_prefix)


if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))
