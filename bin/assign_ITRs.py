#!/usr/bin/env python3
import argparse, sys
from collections import defaultdict

max_is_length = 3000
min_is_length = 700
read_length = 100

def create_cluster_dictionary(itr_clusters):

    seq_num = 0
    with open(itr_clusters, "r") as fa:
        for line in fa:
            if line[0] != ">":
                curr_seq_num = int(line.split(">Seq")[1].split("_")[0])
                if curr_seq_num > seq_num:
                    seq_num = curr_seq_num
    alloc1 = [{}]*seq_num
    alloc2 = [{}]*seq_num

    with open(itr_clusters, "r") as cl:
        for line in cl:
            if line[0] == ">":
                cluster = line.split(" ")[1].split("\n")[0]
            else:
                seq_no = int(line.split(">Seq")[1].split("_")[0]) - 1
                name = line.split("_nstart_")[1].split("_nend_")[0]
                if '_f1' in line:
                    if alloc1[seq_no]:
                        tmp = alloc1[seq_no]
                        tmp[name] = cluster
                        alloc1[seq_no] = tmp
                    else:
                        alloc1[seq_no] = {name: cluster}
                elif '_f2' in line:
                    if alloc2[seq_no]:
                        tmp = alloc2[seq_no]
                        tmp[name] = cluster
                        alloc2[seq_no] = tmp
                    else:
                        alloc2[seq_no] = {name: cluster}
    return((alloc1, alloc2))


def find_itr_clusters(cl_dict, tab_info_file, output_prefix):

    first_contig = ''
    clusters_positions = {}
    itr_clusters = defaultdict(lambda: 0)
    alloc1 = [0]*len(cl_dict[0])
    alloc2 = [0]*len(cl_dict[0])

    with open(output_prefix + '_contigs_reads_itr_position_info.tab', 'w') as out:
        out.write("sample_id\tcontig\tread1\tread2\tposition1\tposition2\titr_cluster\n")
        with open(tab_info_file, "r") as tab:
            next(tab)
            for line in tab:
                seq = line.split('\t')[2]
                current_pos = int(line.split('\t')[3])
                seq_no = int(seq.replace('Seq', '').split('_')[0]) - 1
                name = seq.split("_nstart_")[1].split("_nend_")[0]
                current_contig = line.split('\t')[1]
                if '_f1_LCoord' in seq or '_f2\t' in line:
                    cluster = cl_dict[0][seq_no][name]
                elif '_f2_LCoord' in seq or '_f1\t' in line:
                    cluster = cl_dict[1][seq_no][name]

                if first_contig == '':
                    clusters_positions[cluster] = [line]
                    first_contig = current_contig
                elif first_contig == current_contig:
                    if cluster in clusters_positions:
                        for item in clusters_positions[cluster]:
                            item_fields = item.replace('\n', '').split('\t')
                            pos = int(item_fields[3])
                            if abs(current_pos - pos) <= (max_is_length - read_length*2) and abs(current_pos - pos) >= (min_is_length - read_length*2):
                                line_fields = line.replace('\n', '').split('\t')
                                read1_no = int(item_fields[2].replace('Seq', '').split('_')[0]) - 1
                                read2_no = int(line_fields[2].replace('Seq', '').split('_')[0]) - 1
                                if '_f1_LCoord' in item_fields[2]:
                                    alloc1[read1_no] = 1
                                elif '_f2_LCoord' in item_fields[2]:
                                    alloc2[read1_no] = 1

                                if '_f1_LCoord' in line_fields[2]:
                                    alloc1[read2_no] = 1
                                elif '_f2_LCoord' in line_fields[2]:
                                    alloc2[read2_no] = 1
                                out.write(item_fields[0] + '\t' + current_contig + '\t' + item_fields[2] + '\t' + line_fields[2] + '\t' + item_fields[3] + '\t' + line_fields[3] + '\t' + cluster + '\n')
                    else:
                        clusters_positions[cluster] = [line]
                else:
                    clusters_positions = {cluster: [line]}
                    first_contig = current_contig

    return ((alloc1, alloc2))


def write_fasta(fasta, itr_flags, output_prefix):
    flag = 0
    with open(output_prefix + "_ITRs.fasta", "w") as out:
        with open(fasta, "r") as fa:
            for line in fa:
                if line[0] == ">":
                    if '_f1_LCoord' in line:
                        seq_no = int(line.replace('>Seq', '').split('_')[0]) - 1
                        if itr_flags[0][seq_no]:
                            out.write(line)
                            flag = 1
                    elif '_f2_LCoord' in line:
                        seq_no = int(line.replace('>Seq', '').split('_')[0]) - 1
                        if itr_flags[1][seq_no]:
                            out.write(line)
                            flag = 1
                else:
                    if flag:
                        out.write(line)
                        flag = 0


def get_arguments():
    parser = argparse.ArgumentParser(description='Get the ITR clusters and information on contigs and reads with ITRs.')
    parser.add_argument('--clipped_reads', '-r', dest='clipped_fasta', required=True,
                        help='Input clipped reads FASTA.', type = str)
    parser.add_argument('--cdhit_cluster_file', '-c', dest='cluster_file', required=True,
                        help='Input .clstr file from CD-HIT.', type = str)
    parser.add_argument('--info_tab_file', '-t', dest='tab_file', required=True,
                        help='Input "contigs_reads_ir_position_info.tab" file.', type = str)
    parser.add_argument('--output_prefix', '-o', dest='output_prefix', required=True,
                        help='Prefix of output files.', type = str)
    return parser


def main(args):

    # Create dictionary of IR clusters
    cl_dict = create_cluster_dictionary(args.cluster_file)

    # Find clusters that belong to ITRs
    itr_flags = find_itr_clusters(cl_dict, args.tab_file, args.output_prefix)

    # Write ITR reads
    write_fasta(args.clipped_fasta, itr_flags, args.output_prefix)

if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))
