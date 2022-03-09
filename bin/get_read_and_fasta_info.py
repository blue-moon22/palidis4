#!/usr/bin/env python3
import argparse, sys, os


def get_itr_clusters(is_annot_file):

    contig_clusters = []
    with open(is_annot_file, 'r') as file:
        next(file)
        for line in file:
            contig_clusters.append(line.split('\t')[1] + line.split('\n')[0].split('\t')[6])

    return contig_clusters


def get_read_headers_of_clusters(contigs_itr_clusters, ir_cluster_tab, output_prefix):

    read_headers = []
    with open(f'{output_prefix}_itr_read_positions_clusters.txt', 'w') as out:
        out.write('sample_id\tcontig\tread\titr_start\titr_end\titr_cluster\n')
        with open(ir_cluster_tab, 'r') as file:
            next(file)
            for line in file:
                if line.split('\t')[1] + line.split('\t')[5] in contigs_itr_clusters:
                    out.write('\t'.join(line.split('\t')[0:6]) + '\n')
                    read_headers.append(line.split('\t')[2])

    return read_headers


def write_itr_reads(read_headers, fasta_file, output_prefix):

    flag = 0
    with open(f'{output_prefix}_ITRs.fasta', 'w') as out:
        with open(fasta_file, 'r') as file:
            for line in file:
                if line[0] == '>':
                    if line.split('\n')[0].split('>')[1] in read_headers:
                        out.write(line)
                        flag = 1
                    else:
                        flag = 0
                else:
                    if flag:
                        out.write(line)


def get_arguments():
    parser = argparse.ArgumentParser(description='Get read info for associated ITR clusters.')
    parser.add_argument('--ir_fasta', '-f', dest='fasta_file', required=True,
                        help='IR FASTA file.', type = str)
    parser.add_argument('--ir_cluster_tab', '-t', dest='tab_file', required=True,
                        help='Input "_reads_itr_clusters.txt" file.', type = str)
    parser.add_argument('--is_annotations_tab', '-a', dest='annotations_file', required=True,
                        help='Input IS annotations with ITR clusters.', type = str)
    parser.add_argument('--output_prefix', '-o', dest='output_prefix', required=True,
                    help='Prefix of output files.', type = str)
    return parser


def main(args):

    # Get ITR clusters from IS annotation file
    contigs_itr_clusters = get_itr_clusters(args.annotations_file)

    # Write reads and clusters that are within ITR clusters and return read headers
    read_headers = get_read_headers_of_clusters(contigs_itr_clusters, args.tab_file, args.output_prefix)

    # Write read fasta containing read headers
    write_itr_reads(read_headers, args.fasta_file, args.output_prefix)


if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))
