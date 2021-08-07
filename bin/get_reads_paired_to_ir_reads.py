#!/usr/bin/env python3
import sys, getopt

def get_ir_headers(itr_file):
    print("Getting ITR headers from pal-mem...")
    itr_headers1 = []
    itr_headers2 = []
    with open(itr_file, "r") as f:
        for line in f:
            if line[0] == ">":
                line = line.split("_LCoord")[0]
                if line[-2:] == "f1":
                    itr_headers1.append(line + "\n")
                else:
                    itr_headers2.append(line + "\n")
    return((itr_headers1, itr_headers2))


def identify_ir_reads(itr_headers, fasta1, fasta2):
    print("Allocating flags...")
    numHeaders = 0
    with open(fasta1, "r") as f1:
        for line1 in f1:
            if line1[0] == '>':
                numHeaders += 1
    alloc1 = [1]*int(numHeaders)
    alloc2 = [1]*int(numHeaders)

    print("Getting all non-ITR IDs from FASTA files...")
    # FASTA file 1
    with open(fasta1, "r") as f:
        for line in f:
            if line[0] == ">":
                if line not in itr_headers[0]:
                    alloc1[int(line.split("_")[0].split(">Seq")[1])-1] = 0

    # FASTA file 2
    with open(fasta2, "r") as f:
        for line in f:
            if line[0] == ">":
                if line not in itr_headers[1]:
                    alloc2[int(line.split("_")[0].split(">Seq")[1])-1] = 0
    return([alloc1, alloc2])


def write_reads(fasta1, fasta2, alloc1, alloc2, out_prefix):
    print("Writing paired non-ITR FASTQ 1...")
    numHeaders = 0
    status = 0
    with open(out_prefix + "_non_ITR_1.fasta", "w") as nonitr1:
        with open(fasta1, "r") as f1:
            for line in f1:
                if line[0] == ">":
                    if alloc2[numHeaders] and not alloc1[numHeaders]:
                        status = 1
                        nonitr1.write(line)
                    else:
                        status = 0
                    numHeaders += 1
                else:
                    if status:
                        nonitr1.write(line)

    print("Writing paired non-ITR FASTQ 2...")
    numHeaders = 0
    status = 0
    with open(out_prefix + "_non_ITR_2.fasta", "w") as nonitr2:
        with open(fasta2, "r") as f2:
            for line in f2:
                if line[0] == ">":
                    if alloc1[numHeaders] and not alloc2[numHeaders]:
                        status = 1
                        nonitr2.write(line)
                    else:
                        status = 0
                    numHeaders += 1
                else:
                    if status:
                        nonitr2.write(line)


def main(argv):
    try:
        opts, args = getopt.getopt(argv, "h")
    except getopt.GetoptError:
        print("get_reads_paired_to_ir_reads.py <ITR_file> <fasta_1> <fasta_2> <out_prefix>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print("get_reads_paired_to_ir_reads.py <ITR_file> <fasta_1> <fasta_2> <out_prefix>")
            sys.exit()
    itr_file = sys.argv[1]
    fasta1 = sys.argv[2]
    fasta2 = sys.argv[3]
    output_prefix = sys.argv[4]

    itr_headers = get_ir_headers(itr_file)
    read_bins = identify_ir_reads(itr_headers, fasta1, fasta2)
    write_reads(fasta1, fasta2, read_bins[0], read_bins[1], output_prefix)


if __name__ == "__main__":
    main(sys.argv[1:])
