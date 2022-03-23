#!/usr/bin/env python3
import argparse, sys

from collections import defaultdict


def get_new_clusters(clstr_file):

    new_clusters = {}

    with open(clstr_file, "r") as file:
        for line in file:
            if line[0] == ">":
                cluster = line.split('\n')[0].replace('>Cluster ', '')
            else:
                seq = line.split('>')[1].split('...')[0]
                new_clusters[seq] = cluster

    return new_clusters


def get_itr_reads(itr_reads_tab):

    contig_reads = defaultdict(lambda: [])

    for tab_file in itr_reads_tab:
        with open(tab_file, "r") as file:
            next(file)
            for line in file:
                sample_id = line.split('\t')[0]
                contig = line.split('\t')[1]
                read = line.split('\t')[2]
                start = int(line.split('\t')[3])
                end = int(line.split('\t')[4])
                old_cluster = line.split('\n')[0].split('\t')[5]
                tmp = contig_reads[f'{sample_id}_{contig}_{old_cluster}']
                tmp.append((read, start, end))
                contig_reads[f'{sample_id}_{contig}_{old_cluster}'] = tmp

    return contig_reads


def assign_new_clusters(new_clusters, contig_reads, is_annot_tab, output_prefix):

    with open(f'{output_prefix}_insertion_sequence_catalog.txt', "w") as out:
        out.write("sample_id\tcontig\titr1_start_position\titr1_end_position\titr2_start_position\titr2_end_position\titr_cluster\n")
        for tab_file in is_annot_tab:
            with open(tab_file, "r") as file:
                next(file)
                for line in file:
                    sample_id = line.split('\t')[0]
                    contig = line.split('\t')[1]
                    old_cluster = line.split('\n')[0].split('\t')[6]
                    itr_clusters = []
                    for items in contig_reads[f'{sample_id}_{contig}_{old_cluster}']:
                        if (items[1] >= int(line.split('\t')[2]) or items[1] < 1) and items[2] <= int(line.split('\t')[3]) or items[1] >= int(line.split('\t')[4]) and (items[2] <= int(line.split('\t')[5]) or items[2] > int(contig.split('_')[3])):
                            itr_clusters.append(new_clusters[items[0]])
                    out.write('\t'.join(line.split('\t')[0:6]) + '\t' + ';'.join(sorted(set(itr_clusters))) + '\n')


def get_arguments():
    parser = argparse.ArgumentParser(description='Generate catalog of insertion sequences.')
    parser.add_argument('--clstr_file', '-c', dest='clstr_file', required=True,
                        help='Input CD-HIT cluster file.', type = str)
    parser.add_argument('--itr_reads_tab', '-r', dest='itr_reads_tab', required=True,
                        help='List of "*_itr_read_positions_clusters.txt".', nargs='*')
    parser.add_argument('--is_annot_tab', '-a', dest='is_annot_tab', required=True,
                        help='List of "*_insertion_sequence_annotations.tab".', nargs='*')
    parser.add_argument('--output_prefix', '-o', dest='output_prefix', required=True,
                        help='Prefix of output tab file.', type = str)
    return parser


def main(args):

    # Get reads with their new clusters
    new_clusters = get_new_clusters(args.clstr_file)

    # Get contigs with their associated reads
    contig_reads = get_itr_reads(args.itr_reads_tab)

    # Get IS annotations with new clusters
    assign_new_clusters(new_clusters, contig_reads, args.is_annot_tab, args.output_prefix)


if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))
