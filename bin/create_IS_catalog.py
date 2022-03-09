#!/usr/bin/env python3
import argparse, sys
import pandas as pd

def assign_cluster_bins(clstr_file):

    df = pd.DataFrame(columns=['read', 'new_cluster'], index = [0])
    ind = 0

    with open(clstr_file, "r") as clstr:
        for line in clstr:
            if line[0] == ">":
                cluster = line.split('\n')[0].replace('>Cluster ', '')
                df.loc[ind, 'new_cluster'] = int(cluster)
            else:
                seq = line.split('>')[1].split('...')[0]
                df.loc[ind, 'new_cluster'] = seq


def get_arguments():
    parser = argparse.ArgumentParser(description='Generate catalog of insertion sequences.')
    parser.add_argument('--clstr_file', '-c', dest='clstr_file', required=True,
                        help='Input CD-HIT cluster file.', type = str)
    parser.add_argument('--output_prefix', '-o', dest='output_prefix', required=True,
                        help='Prefix of output tab file.', type = str)
    return parser


def main(args):
    assign_cluster_bins(args.clstr_file)


if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))
