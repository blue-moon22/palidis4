#!/usr/bin/env python3

"""
author: Victoria Carr
Function for clipping reads.
"""

import argparse, sys

def get_positions(tab_file):

    positions = {}
    with open(tab_file, "r") as f:
        for line in f:
            col1 = line.split('\t')[0]
            col2 = line.split('\t')[1]
            header1 = col1.split('_LCoord')[0]
            header2 = col2.split('_LCoord')[0]
            start1 = int(col1.split('_LCoord_')[1].split('_RCoord')[0])
            end1 = int(col1.split('_RCoord_')[1])
            start2 = int(col2.split('_LCoord_')[1].split('_RCoord')[0])
            end2 = int(col2.split('\n')[0].split('_RCoord_')[1])

            if header1 not in positions:
                positions[header1] = []
            tmp = positions[header1]
            tmp.append((start1, end1))
            positions[header1] = tmp

            if header2 not in positions:
                positions[header2] = []
            tmp = positions[header2]
            tmp.append((start2, end2))
            positions[header2] = tmp

    return positions


def clip_reads(fasta_file, positions, output_prefix):
    """
    This function clips reads based on the left-hand and right-hand coordinates
    in the headers
    """

    count = 1
    with open(output_prefix + '_irs.fasta', "w") as out:
        with open(fasta_file, 'r') as f:
            for line in f:
                if line[0] == '>':
                    header = line.replace(">", "").replace("\n", "")
                else:
                    for pos in positions[header]:
                        out.write(f'>{str(count)}\n{line[pos[0]-1:pos[1]]}\n')
                        count += 1


def get_arguments():
    parser = argparse.ArgumentParser(description='Get the inverted repeat sequences from reads containing them.')
    parser.add_argument('--read_fasta', '-f', dest='read_fasta', required=True,
                        help='Input read FASTA file.', type = str)
    parser.add_argument('--tab_file', '-t', dest='tab_file', required=True,
                        help='Input tab file from palmem.', type = str)
    parser.add_argument('--output_prefix', '-o', dest='output_prefix', required=True,
                        help='Prefix of output FASTA file.', type = str)
    return parser


def main(args):

    positions = get_positions(args.tab_file)

    clip_reads(args.read_fasta, positions, args.output_prefix)


if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))
