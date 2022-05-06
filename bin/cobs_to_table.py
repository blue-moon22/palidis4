#!/usr/bin/env python3

import argparse

def get_query_size(fai_file, kmer_length):
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
    with open(outfile, "w") as f:
        f.write("query\tsample_id\tkmer_similarity\n")
        for query, info in results.items():
            for sample, hit in info['hits'].items():
                f.write(query + "\t" + sample + "\t" + str(hit) + "\n")
    return



def init():
    ''' get all the input from the user'''
    parser = argparse.ArgumentParser(
        description='producing a tabular results for the cobs hits',
           )
   # input options
    parser.add_argument('--cobs_outfile',
                        required=True,
                        type=str,
                        help='cobs out file')
    parser.add_argument('--fai_file',
                        required=True,
                        help='samtools index of query file for cobs',
                        type=str)
    parser.add_argument('--outname',
                        required=True,
                        help ='name of table output',
                        type=str)
    parser.add_argument('--kmer_length',  # length of query kmer
                        required=False,
                        type=int, default=31,
                        help='Length of query kmer [%(default)s]')
    options = parser.parse_args()
    return options


if __name__ == "__main__":
    options = init()
    query_info = get_query_size(options.fai_file, options.kmer_length)
    results = get_results(query_info, options.cobs_outfile)
    write_table(results, options.outname)
