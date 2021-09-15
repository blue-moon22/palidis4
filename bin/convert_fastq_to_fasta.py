#!/usr/bin/env python3
import argparse, sys, os
from Bio import SeqIO


def write_fasta(fastq_file, read):
    base_name = os.path.basename(fastq_file).split(".")[0]
    reversed = base_name[::-1]
    find = "_" + str(read)
    sample_name = reversed.replace(find[::-1], "", 1)[::-1]
    fasta_file = sample_name + "_" + str(read) + ".fasta"
    index = 1
    with open(fasta_file, "w") as out:
        for record in SeqIO.parse(fastq_file, "fastq"):
            new_id = ">Seq" + str(index) + "_nstart_" + sample_name + "_nend_" + record.id.replace(" ", "_") + "_f" + str(read)
            out.write(new_id + "\n" + str(record.seq) + "\n")
            index += 1


def get_arguments():
    parser = argparse.ArgumentParser(description='Convert FASTQ to FASTA with appropriate headers.')
    parser.add_argument('--fastq', '-f', dest='fastq', required=True,
                        help='Input FASTQ file.', type = str)
    parser.add_argument('--read', '-r', dest='read', required=True,
                        help='Read file number (1 or 2).', type = str)
    return parser


def main(args):
    write_fasta(args.fastq, args.read)


if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))
