#!/usr/bin/env python3

"""
author: Victoria Carr

Functions to get insertion sequence information.
"""

import argparse, sys
import re
from collections import defaultdict

class irClusterInfo:

    def __init__(self):
        self.ir_clusters = {}

    def load_ir_clusters(self, cluster_file):
        with open(cluster_file, "r") as fa:
            for line in fa:
                if line[0] == ">":
                    cluster = line.split(" ")[1].split("\n")[0]
                else:
                    ir = line.split(">")[1].split("...")[0]
                    self.ir_clusters[ir] = cluster

    def get_ir_cluster(self, ir):
        return(self.ir_clusters[ir])


class contigTransposaseInfo:

    def __init__(self):
        self.contig_info = {}

    def load_contig_info(self, contig_info_file):
        with open(contig_info_file, "r") as file:
            next(file)
            for line in file:
                protein = line.split("\t")[0]
                accession = line.split("\t")[1]
                annotation = line.split("\t")[2]
                start = int(line.split("\t")[3])
                end = int(line.split("\t")[4].replace("\n", ""))
                contig = protein[:protein.rindex('_')]

                if contig not in self.contig_info:
                    self.contig_info[contig] = []
                tmp = self.contig_info[contig]
                tmp.append([accession, annotation, start, end])
                self.contig_info[contig] = tmp

    def get_contig_info(self):
        return(self.contig_info)

    def check_irs_flanking(self, contig, ir_positions, MIN_IS_LEN, MAX_IS_LEN):
        is_info = {}
        for index1, ir_positions1 in enumerate(ir_positions[:-1]):
            for index2, ir_positions2 in enumerate(ir_positions[index1+1:]):
                if ir_positions2[0] > ir_positions1[0]:
                    length = (ir_positions2[1] + 1) - ir_positions1[0]
                    itr1_start = ir_positions1[0]
                    itr1_end = ir_positions1[1]
                    itr2_start = ir_positions2[0]
                    itr2_end = ir_positions2[1]
                else:
                    length = (ir_positions1[1] + 1) - ir_positions2[0]
                    itr1_start = ir_positions2[0]
                    itr1_end = ir_positions2[1]
                    itr2_start = ir_positions1[0]
                    itr2_end = ir_positions1[1]

                if length >= MIN_IS_LEN and length <= MAX_IS_LEN:
                    for transposase in self.contig_info[contig]:
                        start = transposase[2]
                        end = transposase[3]
                        if start > itr1_end and end < itr2_start:
                            if (itr1_start, itr2_end) not in is_info:
                                is_info[(itr1_start, itr2_end)] = []
                            tmp = is_info[(itr1_start, itr2_end)]
                            tmp.append(transposase)
                            is_info[(itr1_start, itr2_end)] = tmp

        return is_info


def process_sam_file(sam_file, ir_cluster_info):
    """
    Function to read information from the sam file and assign positions of IRs
    """

    with open(sam_file, "r") as sam:

        mapping_info = {}
        for line in sam:
            contig = line.split("\t")[2]
            ir = line.split("\t")[0]
            index = int(ir.split("_")[0].replace("Seq", ""))
            pos = int(line.split("\t")[3])
            len_fragment = len(line.split("\t")[9])
            ir_positions = [pos, (pos + len_fragment)-1]

            cluster = ir_cluster_info.get_ir_cluster(ir)

            if contig in mapping_info:
                if cluster not in mapping_info[contig]:
                    mapping_info[contig][cluster] = []
            else:
                mapping_info[contig] = {cluster: []}
            tmp = mapping_info[contig][cluster]
            tmp.append(ir_positions)
            mapping_info[contig][cluster] = tmp

    return mapping_info


def write_insertion_sequences_info(mapping_info, contig_transposase_info, sample_id, MIN_IS_LEN, MAX_IS_LEN):

    contig_is_info = {}
    with open(f'{sample_id}_insertion_sequences_info.txt', "w") as out:
        out.write("IS_name\tsample_id\tcontig\tis_start\tis_end\tdescription\n")
        for contig, clusters in mapping_info.items():
            for cluster, itr_positions in clusters.items():
                is_info = contig_transposase_info.check_irs_flanking(contig, itr_positions, MIN_IS_LEN, MAX_IS_LEN)
                for is_position, transposases in is_info.items():
                    length = (is_position[1]- is_position[0]) + 1
                    is_name = f'IS_length_{length}'
                    description = ''
                    for transposase in transposases:
                        is_name = (f'{is_name}_{transposase[0]}_{transposase[2]}_{transposase[3]}')
                        if description:
                            description = (f'{description};{transposase[1]}')
                        else:
                            description = (f'{description}{transposase[1]}')
                    out.write(f'{is_name}\t{sample_id}\t{contig}\t{is_position[0]}\t{is_position[1]}\t{description}\n')
                    if contig not in contig_is_info:
                        contig_is_info[contig] = []
                    tmp = contig_is_info[contig]
                    tmp.append([is_name, is_position[0], is_position[1]])
                    contig_is_info[contig] = tmp

    return contig_is_info


def write_insertion_sequences_fasta(contig_is_info, contig_fasta, output_prefix):

    contig = None
    with open(f'{output_prefix}_insertion_sequences.fasta', "w") as out:
        with open(contig_fasta, "r") as fasta:
            for line in fasta:
                if line[0] == ">":
                    contig = line.replace(">", "").replace("\n", "")
                else:
                    if contig in contig_is_info:
                        for info in contig_is_info[contig]:
                            out.write(f'>{info[0]}\n{line[info[1]-1:(info[2])]}\n')


def get_arguments():

    parser = argparse.ArgumentParser(description='Get reads that contain candidate ITRs and contigs that associate with these.')
    parser.add_argument('--contig_fasta', '-f', dest='contig_fasta', required=True,
                        help='Input contig FASTA file.', type = str)
    parser.add_argument('--contig_info', '-i', dest='contig_info', required=True,
                        help='Input txt file with transposase information.', type = str)
    parser.add_argument('--sam_file', '-s', dest='sam_file', required=True,
                        help='Input sam file.', type = str)
    parser.add_argument('--clstr_file', '-c', dest='clstr_file', required=True,
                        help='Input cluster file.', type = str)
    parser.add_argument('--min_is_len', '-min', dest='min_is_len', required=False,
                        help='Minimum length of insertion sequence', type = int, default = 500)
    parser.add_argument('--max_is_len', '-max', dest='max_is_len', required=False,
                        help='Maximum length of insertion sequence', type = int, default = 3000)
    parser.add_argument('--output_prefix', '-o', dest='output_prefix', required=True,
                        help='Prefix of output files (should be sample ID).', type = str)
    return parser


def main(args):

    # Load cluster information
    ir_cluster_info = irClusterInfo()
    ir_cluster_info.load_ir_clusters(args.clstr_file)

    # Get IR mapping info
    mapping_info = process_sam_file(args.sam_file, ir_cluster_info)

    # Get transposase info
    contig_transposase_info = contigTransposaseInfo()
    contig_transposase_info.load_contig_info(args.contig_info)

    # Write insertion sequence info
    contig_is_info = write_insertion_sequences_info(mapping_info, contig_transposase_info, args.output_prefix, args.min_is_len, args.max_is_len)

    # Write insertion sequence fasta
    write_insertion_sequences_fasta(contig_is_info, args.contig_fasta, args.output_prefix)


if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))
