#!/usr/bin/env python3

"""
author: Victoria Carr and Grace Blackwell
email: victoria.carr@sanger.ac.uk

Functions for tabulating the output of COBS query.
"""

import argparse, sys

def get_query_size(fai_file, kmer_length):
    """
    This function returns a dictionary the maximum number of kmers for each query
    """

    sizes = {}
    with open(fai_file, "r") as f:
        for line in f:
            l = line.strip().split("\t")
            query = l[0]
            size = int(l[1])
            kmers = size - kmer_length + 1
            if query not in sizes:
                sizes[query] = {}
                sizes[query]['max_kmers'] = kmers
                sizes[query]['hits'] = {}
    return sizes


def get_results(queries, results):
    """
    This function returns a dictionary with the alignment identity for each
    query and hit
    """

    with open(results, "r") as f:
        for line in f:
            if line.startswith("*"):
                query = line.strip().split("\t")[0].split("*")[1]
                #check how they are separated and how that matches up with the .fai name
                q = query.split(" ")[0]
            else:
               l = line.strip().split("\t")
               sample = l[0]
               kmers = int(l[1])
               if sample not in queries[q]['hits']:
                   queries[q]['hits'][sample] = kmers/queries[q]['max_kmers']
    return queries


def write_table(results, outfile):
    """
    This function writes the dictionary of information as tab-delimited file
    """

    with open(outfile, "w") as f:
        f.write("query\tsample_id\tkmer_similarity\n")
        for query, info in results.items():
            for sample, hit in info['hits'].items():
                f.write(query + "\t" + sample + "\t" + str(hit) + "\n")


def get_arguments():
    parser = argparse.ArgumentParser(description='Producing a tabular results for the cobs hits.')
    parser.add_argument('--cobs_outfile', required=True, type=str,
                        help='cobs out file')
    parser.add_argument('--fai_file', required=True,
                        help='samtools index of query file for cobs', type=str)
    parser.add_argument('--outname', required=True, type=str,
                        help ='name of table output')
    parser.add_argument('--kmer_length', required=False, type=int, default=31,
                        help='Length of query kmer [%(default)s]')
    return parser


def main(args):

    query_info = get_query_size(args.fai_file, args.kmer_length)

    results = get_results(query_info, args.cobs_outfile)

    write_table(results, args.outname)


if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))
