#!/usr/bin/env python3
import argparse, sys

def create_tab_file(clstr_file, output_prefix):
    with open(output_prefix + ".tab", "w") as out:
        out.write("seq\tcd_hit_cluster\n")
        with open(clstr_file, "r") as clstr:
            for line in clstr:
                if line[0] == ">":
                    cluster = line.split('\n')[0].replace('>Cluster ', '')
                else:
                    seq = line.split('>')[1].split('...')[0]
                    out.write(seq + '\t' + cluster + '\n')


def get_arguments():
    parser = argparse.ArgumentParser(description='Get the inverted repeat sequences from reads containing them.')
    parser.add_argument('--clstr_file', '-c', dest='clstr_file', required=True,
                        help='Input CD-HIT cluster file.', type = str)
    parser.add_argument('--output_prefix', '-o', dest='output_prefix', required=True,
                        help='Prefix of output tab file.', type = str)
    return parser


def main(args):
    create_tab_file(args.clstr_file, args.output_prefix)


if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))
