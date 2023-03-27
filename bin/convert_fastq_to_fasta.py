#!/usr/bin/env python3

"""
author: Victoria Carr

Function to convert a FASTQ file into a FASTA file.
"""

import argparse, sys, os

def write_fasta(fastq_file, read, fasta_file):
    """
    This function writes a FASTA file from a FASTQ file
    """
    index = 1
    count = 1
    with open(fasta_file, "w") as out:
        with open(fastq_file, "r") as fastq:
            for line in fastq:
                if count%4 == 1:
                    new_id = ">Seq" + str(index) + "_f" + str(read)
                    out.write(new_id + "\n")
                    index += 1
                elif count%4 == 2:
                    out.write(line)
                count += 1


def get_arguments():
    parser = argparse.ArgumentParser(description='Convert FASTQ to FASTA with appropriate headers.')
    parser.add_argument('--fastq', '-f', dest='fastq', required=True,
                        help='Input FASTQ file.', type = str)
    parser.add_argument('--fasta', '-o', dest='fasta', required=True,
                        help='Output FASTA file.', type = str)
    parser.add_argument('--read', '-r', dest='read', required=True,
                        help='Read file number (1 or 2).', type = str)
    return parser


def main(args):
    write_fasta(args.fastq, args.read, args.fasta)


if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))
