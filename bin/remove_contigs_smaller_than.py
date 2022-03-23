#!/usr/bin/env python3
import argparse, sys

def select_contigs(fasta_file, min_length, output_file):
    seq = ''
    with open(output_file, 'w') as out:
        with open(fasta_file, 'r') as f:
            for line in f:
                if line[0] == '>':
                    if len(seq) >= min_length:
                        out.write(header + seq + '\n')
                    header = line
                    seq = ''
                else:
                    seq += line.replace('\n', '')


def get_arguments():
    parser = argparse.ArgumentParser(description='Remove contigs less than specified.')
    parser.add_argument('--fasta_file', '-i', dest='fasta_file', required=True,
                        help='Input assembly/contig FASTA file.', type = str)
    parser.add_argument('--output_file', '-o', dest='output_file', required=True,
                        help='Output FASTA file.', type = str)
    parser.add_argument('--min_length', '-l', dest='min_length', required=True,
                        help='Input length.', type = int)
    return parser


def main(args):
    select_contigs(args.fasta_file, args.min_length, args.output_file)


if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))
