#!/usr/bin/env python3
import argparse, sys


def get_read_names(read_names_file):
    read_names = []
    with open(read_names_file, "r") as txt:
        for line in txt:
            read_names.append(">" + line)
    return read_names


def write_fasta(read_names, fasta, output_fasta):
    flag = 0
    with open(output_fasta, "w") as out:
        with open(fasta, "r") as fa:
            for line in fa:
                if line[0] == ">":
                    if line in read_names:
                        out.write(line)
                        flag = 1
                else:
                    if flag:
                        out.write(line)
                        flag = 0


def get_arguments():
    parser = argparse.ArgumentParser(description='Get the ITR clusters and information on contigs and reads with ITRs.')
    parser.add_argument('--input_fasta', '-i', dest='fasta', required=True,
                        help='Input FASTA file.', type = str)
    parser.add_argument('--read_names', '-r', dest='read_names_file', required=True,
                        help='Text file of read names.', type = str)
    parser.add_argument('--output_fasta', '-o', dest='output_fasta', required=True,
                        help='Output FASTA file name.', type = str)
    return parser


def main(args):

    # Get read names
    read_names = get_read_names(args.read_names_file)

    # Write output FASTA file
    write_fasta(read_names, args.fasta, args.output_fasta)


if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))
