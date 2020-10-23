#!/usr/bin/env python3

import sys, getopt

def write_fasta(fasta1, fasta2, itr_file, out_prefix):

    print("Allocating flags...")
    numHeaders = 0
    with open(fasta1, "r") as f1:
        for line1 in f1:
            numHeaders += 1
    alloc1 = [0]*int(numHeaders/2)
    alloc2 = [0]*int(numHeaders/2)

    print("Reading ITR file...")
    with open(itr_file, "r") as itr:
        for line in itr:
            if line[0] == ">":
                line = line.split("_LCoord")[0]
                if line[-2:] == "f1":
                    alloc2[int(line.split("_")[0].replace(">Seq", ""))-1] = 1
                else:
                    alloc1[int(line.split("_")[0].replace(">Seq", ""))-1] = 1

    print("Writing discordant reads...")
    flag = 0
    with open(out_prefix + "_discord_1.fasta", "w") as d1:
        with open(fasta1, "r") as f1:
            for line in f1:
                if line[0] == ">":
                    if alloc1[int(line.split("_")[0].replace(">Seq", ""))-1] and not alloc2[int(line.split("_")[0].replace(">Seq", ""))-1]:
                        d1.write(line)
                        flag = 1
                    else:
                        flag = 0
                else:
                    if flag:
                        d1.write(line)

    flag = 0
    with open(out_prefix + "_discord_2.fasta", "w") as d2:
        with open(fasta2, "r") as f2:
            for line in f2:
                if line[0] == ">":
                    if alloc2[int(line.split("_")[0].replace(">Seq", ""))-1] and not alloc1[int(line.split("_")[0].replace(">Seq", ""))-1]:
                        d2.write(line)
                        flag = 1
                    else:
                        flag = 0
                else:
                    if flag:
                        d2.write(line)

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "h")
    except getopt.GetoptError:
        print("get_discordant_reads.py <fasta_file_1> <fasta_file_2> <ITR_file> <out_prefix>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print("get_discordant_reads.py <fasta_file_1> <fasta_file_2> <ITR_file> <out_prefix>")
            sys.exit()
    fasta1 = sys.argv[1]
    fasta2 = sys.argv[2]
    itr_file = sys.argv[3]
    out_prefix = sys.argv[4]

    write_fasta(fasta1, fasta2, itr_file, out_prefix)


if __name__ == "__main__":
    main(sys.argv[1:])
