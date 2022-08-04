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


def write_info(tab_file, cobs_info, output_prefix):
    """
    Function to write the annotation information of the insertion sequences
    """

    with open(f'{output_prefix}_insertion_sequences_info.txt', "w") as out:
        out.write("IS_name\tsample_id\tcontig\titr1_start_position\titr1_end_position\titr2_start_position\titr2_end_position\titr_cluster\tCOBS_index_biosample_id\tCOBS_index_origin\n")
        with open(tab_file, "r") as file:
            next(file)
            for line in file:
                out.write(line.replace('\n', ''))
                is_name = line.split('\t')[0]
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

    parser = argparse.ArgumentParser(description='Get tab file of annotations.')
    parser.add_argument('--tab_file', '-t', dest='tab_file', required=True,
                        help='Input "insertion_sequence_annotations.tab" file.', type = str)
    parser.add_argument('--cobs_search_out', '-c', dest='cobs_table', required=False,
                        help='Input "_results_table.txt" file.', type = str)
    parser.add_argument('--ffq_json', '-j', dest='ffq_json', required=False,
                        help='Location of JSON files.', type = str)
    parser.add_argument('--output_prefix', '-o', dest='output_prefix', required=True,
                    help='Prefix of output files.', type = str)
    return parser


def main(args):

    if args.cobs_table and args.ffq_json:
        cobs_info_dict = get_cobs_info(args.cobs_table, args.ffq_json)
    else:
        cobs_info_dict = {}

    write_info(args.tab_file, cobs_info_dict, args.output_prefix)


if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))
