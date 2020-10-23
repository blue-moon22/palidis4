#!/usr/bin/env python3

import sys, getopt

def write_fasta(fasta, sam_file, out_fasta):

    contigs = []
    with open(sam_file, "r") as sam:
        for line in sam:
            contigs.append(line.split("\t")[2])

    flag = 0
    with open(out_fasta, "w") as out:
        with open(fasta, "r") as f:
            for line in f:
                if line[0] == ">":
                    if line.split("\n")[0].replace(">", "") in contigs:
                        out.write(line)
                        flag = 1
                    else:
                        flag = 0
                else:
                    if flag:
                        out.write(line)

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "h")
    except getopt.GetoptError:
        print("get_contigs_with_sites.py <contig_fasta> <contig_mapping_sam_file> <out_fasta>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print("get_contigs_with_sites.py <contig_fasta> <contig_mapping_sam_file> <out_fasta>")
            sys.exit()
    fasta = sys.argv[1]
    sam_file = sys.argv[2]
    out_fasta = sys.argv[3]

    write_fasta(fasta, sam_file, out_fasta)


if __name__ == "__main__":
    main(sys.argv[1:])
