#!/usr/bin/env python3

"""
author: Victoria Carr
email: victoria.carr@sanger.ac.uk

Function for clipping reads.
"""

import argparse, sys

def clip_reads(fasta_file, output_prefix):
    """
    This function clips reads based on the left-hand and right-hand coordinates
    in the headers
    """

    with open(output_prefix + '_irs.fasta', "w") as out:
        with open(fasta_file, 'r') as f:
            for line in f:
                if line[0] == '>':
                    out.write(line)
                    start = int(line.split('_LCoord_')[1].split('_RCoord')[0]) - 1
                    end = int(line.split('\n')[0].split('_RCoord_')[1])
                else:
                    out.write(line[start:end] + '\n')


def get_arguments():
    parser = argparse.ArgumentParser(description='Get the inverted repeat sequences from reads containing them.')
    parser.add_argument('--read_fasta', '-f', dest='read_fasta', required=True,
                        help='Input read FASTA file.', type = str)
    parser.add_argument('--output_prefix', '-o', dest='output_prefix', required=True,
                        help='Prefix of output FASTA file.', type = str)
    return parser


def main(args):
    clip_reads(args.read_fasta, args.output_prefix)


if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))
