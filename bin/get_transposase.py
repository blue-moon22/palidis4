#!/usr/bin/env python3

import sys, getopt
import pickle

def write_fasta(hmm_out, sam_file, fasta, out_prefix):

    # Get transposase contigs and positions
    trsps_line = []
    trsps_contig = []
    trsps_positions = []
    with open(hmm_out, "r") as b:
        for line in b:
            if line[:4] == "NODE":
                trsps_line.append(line)
                contig = line.split(" ")[0]
                trsps_contig.append("_".join(contig.split("_")[:len(contig.split("_"))-1]))
                pos1 = int(line.split("#")[1].split(" ")[1])
                pos2 = int(line.split("#")[2].split(" ")[1])
                if pos1 > pos2:
                    trsps_positions.append((pos2, pos1))
                else:
                    trsps_positions.append((pos1, pos2))

    # Write contig tab file
    contigs_with_is = []
    with open(out_prefix + "_contig_is_sites.tab", "w") as out:
        out.write("contig" + "\t" + "transposase_seq" + "\t" + "transposase_start" + "\t" + "transposase_end" + "\t" + "itr_position" + "\n")
        with open(sam_file, "r") as sam:
            for line in sam:
                node = line.split("\t")[2]
                pos = int(line.split("\t")[3])
                for ind, contig in enumerate(trsps_contig):
                    if node == contig:
                        if abs(pos - trsps_positions[ind][0]) < 2750 and abs(pos - trsps_positions[ind][1]) < 2750:
                            contigs_with_is.append(contig)
                            out.write(contig + "\t" + trsps_line[ind].split(" ")[0] + "\t" + str(trsps_positions[ind][0]) + "\t" + str(trsps_positions[ind][1]) + "\t" + str(pos) + "\n")

    # Write contig fasta file
    flag = 0
    with open(out_prefix + "_contig_is_sites.fasta", "w") as out:
        with open(fasta, "r") as f:
            for line in f:
                if line[0] == ">":
                    if line.split("\n")[0].replace(">", "") in contigs_with_is:
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
        print("get_transposase.py <transposase_hmm_out> <contig_mapping_sam_file> <contig_fasta> <out_prefix>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print("get_transposase.py <transposase_hmm_out> <contig_mapping_sam_file> <contig_fasta> <out_prefix>")
            sys.exit()
    hmm_out = sys.argv[1]
    sam_file = sys.argv[2]
    fasta = sys.argv[3]
    out_prefix = sys.argv[4]

    write_fasta(hmm_out, sam_file, fasta, out_prefix)


if __name__ == "__main__":
    main(sys.argv[1:])
