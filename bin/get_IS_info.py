#!/usr/bin/env python3

"""
author: Victoria Carr
email: victoria.carr@sanger.ac.uk

Functions to get the information for insertion sequences.
"""

import argparse, sys, os
import json
import re

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


def write_info(tab_file, prodigal_info, interpro_info, cobs_info, output_prefix):
    """
    Function to write the annotation information of the insertion sequences
    """

    with open(f'{output_prefix}_insertion_sequences_info.txt', "w") as out:
        out.write("IS_name\tsample_id\tcontig\titr1_start_position\titr1_end_position\titr2_start_position\titr2_end_position\titr_cluster\tinterpro_or_panther_accession\tCOBS_index_biosample_id\tCOBS_index_origin\n")
        with open(tab_file, "r") as file:
            next(file)
            for line in file:
                is_name = line.split('\t')[0]

                # Interpro info
                if is_name in prodigal_info:
                    flag = 0
                    proteins = prodigal_info[is_name].keys()
                    for protein in proteins:
                        if protein in interpro_info:
                            accessions = []
                            for acc, items in interpro_info[protein].items():
                                length = is_name.split('_')[4]
                                transposase = '_'.join(re.split('-| ', items[0]))
                                transposase = ''.join(re.split('[^a-zA-Z0-9_]*', transposase))
                                protein_start = int(prodigal_info[is_name][protein][0])
                                start = str((protein_start - 1) + items[1][0])
                                end = str((protein_start - 1) + items[1][1])
                                if flag:
                                    out.write(f'_{transposase}_{start}-{end}')
                                else:
                                    out.write(f'IS_length_{length}_{transposase}_{start}-{end}')
                                    flag = 1
                                accessions.append(acc)

                            out.write('\t')

                    if flag:
                        out.write('\t'.join(line.replace('\n', '').split('\t')[1:]) + '\t' + ';'.join(accessions))

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


def get_prodigal_info(aa_fasta):

    prodigal_dict = {}
    with open(aa_fasta) as faa:
        for line in faa:
            if line[0] == '>':
                header = line.split(' ')[0].replace('>', '')
                contig = header.rsplit('_', 1)[0]
                start = int(line.split(' ')[2])
                end = int(line.split(' ')[4])
                if contig in prodigal_dict:
                    prodigal_dict[contig][header] = (start, end)
                else:
                    prodigal_dict[contig] = {header: (start, end)}

    return prodigal_dict


def get_interproscan_info(interproscan_out):

    interpro_dict = {}
    with open(interproscan_out) as tsv:
        for line in tsv:
            if line.split('\t')[3] == 'PANTHER':
                annotation = line.split('\t')[5]
            else:
                annotation = line.split('\t')[12].replace('\n', '')

            if any(x in ''.join(re.split('[^a-zA-Z0-9]*', annotation)).lower() for x in ['transposase', 'integraselike', 'ribonucleaseh']):
                protein = line.split('\t')[0]
                start = (int(line.split('\t')[6])-1)*3
                end = (int(line.split('\t')[7])-1)*3

                if line.split('\t')[3] == 'PANTHER':
                    accession = line.split('\t')[4]
                else:
                    accession = line.split('\t')[11]

                if protein in interpro_dict:
                    interpro_dict[protein][accession] = [annotation, (start, end)]
                else:
                    interpro_dict[protein] = {accession: [annotation, (start, end)]}

    return interpro_dict


def get_arguments():

    parser = argparse.ArgumentParser(description='Get tab file of annotations.')
    parser.add_argument('--tab_file', '-t', dest='tab_file', required=True,
                        help='Input "insertion_sequence_annotations.tab" file.', type = str)
    parser.add_argument('--cobs_search_out', '-c', dest='cobs_table', required=False,
                        help='Input "_results_table.txt" file.', type = str)
    parser.add_argument('--ffq_json', '-j', dest='ffq_json', required=False,
                        help='Location of JSON files.', type = str)
    parser.add_argument('--aa_fasta', '-p', dest='aa_fasta', required=True,
                        help='Input "_insertion_sequences.faa file."', type = str)
    parser.add_argument('--interproscan_out', '-i', dest='interproscan_out', required=True,
                        help='Input "_insertion_sequences.faa.tsv" file."', type = str)
    parser.add_argument('--output_prefix', '-o', dest='output_prefix', required=True,
                    help='Prefix of output files.', type = str)
    return parser


def main(args):

    prodigal_info_dict = get_prodigal_info(args.aa_fasta)
    interpro_info_dict = get_interproscan_info(args.interproscan_out)

    if args.cobs_table and args.ffq_json:
        cobs_info_dict = get_cobs_info(args.cobs_table, args.ffq_json)
    else:
        cobs_info_dict = {}

    write_info(args.tab_file, prodigal_info_dict, interpro_info_dict, cobs_info_dict, args.output_prefix)


if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))
