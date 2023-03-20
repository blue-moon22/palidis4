#!/usr/bin/env python3

"""
author: Victoria Carr

Function to output contigs with transposase.
"""

import argparse, sys
import re

def get_interproscan_info(interproscan_info, output_prefix):

    contigs = []
    with open(f'{output_prefix}_filtered_transposase.tsv', "w") as out:
        out.write("protein_name\taccession\tannotation\tnucl_start\tnucl_end\n")
        with open(interproscan_info) as tsv:
            for line in tsv:

                panther = 0
                annotation = line.split('\t')[12].replace('\n', '')

                if annotation == '-' and line.split('\t')[3] == 'PANTHER':
                    annotation = line.split('\t')[5]
                    panther = 1

                if any(x in ''.join(re.split('[^a-zA-Z0-9]*', annotation)).lower() for x in ['transposase', 'integraselike', 'ribonucleaseh']):
                    protein = line.split('\t')[0]
                    start = (int(line.split('\t')[6]))*3-2
                    end = (int(line.split('\t')[7]))*3

                    if panther:
                        accession = line.split('\t')[4]
                    else:
                        accession = line.split('\t')[11]

                    contig = protein[:protein.rindex('_')]
                    if contig not in contigs:
                        contigs.append(contig)
                    out.write(f'{protein}\t{accession}\t{annotation}\t{start}\t{end}\n')

    return contigs

def get_contigs(fasta_file, contigs, output_prefix):

    flag = 0
    with open(f'{output_prefix}_filtered_transposase.fasta', "w") as out:
        with open(fasta_file) as fasta:
            for line in fasta:
                if line[0] == ">":
                    contig = line.replace("\n", "").replace(">", "")
                    if contig in contigs:
                        flag = 1
                        out.write(line)
                    else:
                        flag = 0
                else:
                    if flag:
                        out.write(line)

def get_arguments():

    parser = argparse.ArgumentParser(description='Get contigs with transposase.')
    parser.add_argument('--fasta_file', '-f', dest='fasta_file', required=True,
                        help='Input assembly/contig FASTA file.', type = str)
    parser.add_argument('--output_prefix', '-o', dest='output_prefix', required=True,
                        help='Output prefix.', type = str)
    parser.add_argument('--interproscan_tsv', '-t', dest='tsv_file', required=True,
                        help='Interproscan tsv.', type = str)
    return parser


def main(args):

    contigs = get_interproscan_info(args.tsv_file, args.output_prefix)
    get_contigs(contigs, args.output_prefix)

if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))
