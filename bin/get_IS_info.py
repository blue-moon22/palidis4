#!/usr/bin/env python3

"""
author: Victoria Carr
email: victoria.carr@sanger.ac.uk

Functions to get the information for insertion sequences.
"""

import argparse, sys, os
import json

ALIGN_THRESHOLD = 0.99
IDENTITY_THRESHOLD = 99

def create_annotations_dict(blast_out):
    """
    Function to return a dictionary containing annotation information
    """

    annot_dict = {}
    with open(blast_out, "r") as file:
        for line in file:
            query_seq = line.split('\t')[0]
            al_len = int(line.split('\t')[3])
            identity = float(line.split('\t')[2])
            len_seq = int(query_seq.split('_')[4])
            subject_seq = line.split('\t')[1]

            if al_len/len_seq >= ALIGN_THRESHOLD and identity >= IDENTITY_THRESHOLD:
                annot = ('IS', subject_seq.split('_')[0])
            else:
                annot = ('IS_family', subject_seq.split('_')[1])

            if query_seq not in annot_dict:
                annot_dict[query_seq] = [annot]
            else:
                if annot_dict[query_seq][0][0] != 'IS':
                    if annot[0] == 'IS':
                        annot_dict[query_seq] = annot
                    else:
                        tmp = annot_dict[query_seq]
                        if annot not in tmp:
                            tmp.append(annot)
                        annot_dict[query_seq] = tmp
                else:
                    if annot[0] == 'IS' and annot not in annot_dict[query_seq]:
                        tmp = annot_dict[query_seq]
                        if annot not in tmp:
                            tmp.append(annot)
                        annot_dict[query_seq] = tmp

    return annot_dict


def get_isfinder_info(is_info_csv):
    """
    Function to return a dictionary with the origin species for insertion
    sequences from ISfinder
    """

    is_info_dict = {}
    with open(is_info_csv, "r") as file:
        next(file)
        for line in file:
            is_name = line.split(',')[1]
            origin = line.split(',')[6]
            is_info_dict[is_name] = origin

    return is_info_dict


def get_cobs_info(cobs_table, json_loc):
    """
    Function to return a dictionary containing the taxonomy of the origin for
    every IS found in the COBS query
    """

    ffq_info = {}

    for file in os.listdir(json_loc):
        if file.endswith('json'):
            f = open(f'{json_loc}/{file}')
            try:
                json_content = json.load(f)
            except:
                json_content = {}
            ffq_info.update(json_content)
            f.close()

    cobs_info_dict = {}

    with open(cobs_table, "r") as file:
        next(file)
        for line in file:
            query = line.split('\t')[0]
            biosample_id = line.split('\t')[1]
            if biosample_id[:4] == 'SAMN' and biosample_id in ffq_info:
                organism = ffq_info[biosample_id]['samples']['organism']
            else:
                organism = "unknown"

            if query not in cobs_info_dict:
                cobs_info_dict[query] = []
            tmp = cobs_info_dict[query]
            tmp.append((biosample_id, organism))
            cobs_info_dict[query] = tmp

    return cobs_info_dict


def write_info(tab_file, is_finder_annot, is_finder_info, cobs_info, output_prefix):
    """
    Function to write the annotation information of the insertion sequences
    """

    with open(f'{output_prefix}_insertion_sequences_info.txt', "w") as out:
        out.write("IS_name\tsample_id\tcontig\titr1_start_position\titr1_end_position\titr2_start_position\titr2_end_position\titr_cluster\tISfinder_name\tISfinder_origin\tpredicted_IS_family\tCOBS_index_biosample_id\tCOBS_index_origin\n")
        with open(tab_file, "r") as file:
            next(file)
            for line in file:
                is_name = line.split('\t')[0]
                # Get ISfinder info
                if is_name in is_finder_annot:
                    is_finder_names = []
                    is_finder_origins = []
                    if is_finder_annot[is_name][0][0] == 'IS':
                        for item in is_finder_annot[is_name]:
                            is_finder_names.append(item[1])
                            is_finder_origins.append(is_finder_info[item[1]])
                    else:
                        for item in is_finder_annot[is_name]:
                            is_finder_names.append(item[1])

                    if is_finder_annot[is_name][0][0] == 'IS':
                        out.write(line.replace('\n', '') + '\t' + ';'.join(is_finder_names) + '\t' + ';'.join(is_finder_origins) + '\t')
                    else:
                        out.write(line.replace('\n', '') + '\t\t\t' + ';'.join(is_finder_names))
                else:
                    out.write(line.replace('\n', '') + '\t\t\t')

                # Get COBS index info
                if is_name in cobs_info:
                    cobs_biosample = []
                    cobs_origin = []
                    for item in cobs_info[is_name]:
                        cobs_biosample.append(item[0])
                        cobs_origin.append(item[1])
                    out.write('\t' + ';'.join(cobs_biosample) + '\t' + ';'.join(cobs_origin) + '\n')
                else:
                    out.write('\t\t\n')


def get_arguments():

    parser = argparse.ArgumentParser(description='Get tab file of ISfinder annotations.')
    parser.add_argument('--blast_out', '-b', dest='blast_out', required=True,
                        help='BLAST output file.', type = str)
    parser.add_argument('--tab_file', '-t', dest='tab_file', required=True,
                        help='Input "insertion_sequence_annotations.tab" file.', type = str)
    parser.add_argument('--is_finder_info', '-i', dest='is_info_csv', required=True,
                        help='Input ".csv" file.', type = str)
    parser.add_argument('--cobs_search_out', '-c', dest='cobs_table', required=False,
                        help='Input "_results_table.txt" file.', type = str)
    parser.add_argument('--ffq_json', '-j', dest='ffq_json', required=False,
                        help='Location of JSON files.', type = str)
    parser.add_argument('--output_prefix', '-o', dest='output_prefix', required=True,
                    help='Prefix of output files.', type = str)
    return parser


def main(args):

    isfinder_annot_dict = create_annotations_dict(args.blast_out)
    isfinder_info_dict = get_isfinder_info(args.is_info_csv)

    if args.cobs_table and args.ffq_json:
        cobs_info_dict = get_cobs_info(args.cobs_table, args.ffq_json)
    else:
        cobs_info_dict = {}

    write_info(args.tab_file, isfinder_annot_dict, isfinder_info_dict, cobs_info_dict, args.output_prefix)


if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))
