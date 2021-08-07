#!/usr/bin/env python3
import argparse, sys
from collections import defaultdict

def get_largest_index(fasta_file):
    """
    Get the largest sequence ID
    """
    with open(fasta_file, "r") as f:
        for line in f:
            if line[0] == '>':
                index = int(line.split("_")[0].replace(">Seq", ""))
    return(index)


def get_array_size(fasta_file1, fasta_file2):
    size = max(get_largest_index(fasta_file1), get_largest_index(fasta_file2))
    return(size)


def allocate_seqs(size):
    """
    Get the offset length of the IRs and assign to sequence ID
    """
    alloc1 = [0]*size
    alloc2 = [0]*size
    return([alloc1, alloc2])


def get_ir_offset(fasta_file1, fasta_file2, size):
    """
    Get the offset length of the IRs and assign to sequence ID
    """
    alloc1 = [0]*size
    with open(fasta_file1, "r") as f1:
        for line in f1:
            if line[0] == ">":
                index = int(line.split("_")[0].replace(">Seq", ""))
                alloc1[index-1] = int(line.split("\t")[0].split('_RCoord')[0].split('LCoord_')[1])

    alloc2 = [0]*size
    with open(fasta_file2, "r") as f2:
        for line in f2:
            if line[0] == ">":
                index = int(line.split("_")[0].replace(">Seq", ""))
                alloc2[index-1] = int(line.split("\t")[0].split('_RCoord')[0].split('LCoord_')[1])
    return([alloc1, alloc2])

def write_tab_file(itr_info1, itr_info2, out_prefix):

    # Write tab file
    with open(out_prefix + '_contigs_reads_ir_position_info.tab', "w") as out_tab:
        out_tab.write('sample_id\tcontig\tread\tposition\n')
        sample_id = out_prefix.split('/')[len(out_prefix.split('/'))-1]
        for contig in itr_info1:
            for read in itr_info1[contig]:
                out_tab.write(sample_id + '\t' + contig + '\t' + read + '\t' + str(itr_info1[contig][read]) + '\n')

        for contig in itr_info2:
            for read in itr_info2[contig]:
                out_tab.write(sample_id + '\t' + contig + '\t' + read + '\t' + str(itr_info2[contig][read]) + '\n')


def write_contig_file(contigs_fasta_file, itr_info1, itr_info2, out_prefix):
    # Write fasta file
    with open(out_prefix + '_contigs_with_candidate_itrs.fa', "w") as out_fasta:
        flag = 0
        with open(contigs_fasta_file, "r") as f:
            for line in f:
                if line[0] == ">":
                    contig = line.split("\n")[0].split(" ")[0].replace(">", "")
                    if contig in itr_info1 or contig in itr_info2:
                        out_fasta.write(line)
                        flag = 1
                    else:
                        flag = 0
                else:
                    if flag:
                        out_fasta.write(line)


def write_fasta_files(read_fasta1, read_fasta2, alloc_arrays, out_prefix):
    # Get read names of IRs from first and second paired files
    # Write fasta file
    alloc1 = alloc_arrays[0]
    alloc2 = alloc_arrays[1]
    read_pairs_dict = {}
    with open(out_prefix + '_reads_with_candidate_itrs_1.fasta', "w") as out_fasta1:
        flag = 0
        with open(read_fasta1, "r") as f1:
            for line in f1:
                if line[0] == ">":
                    index = int(line.split("_")[0].replace(">Seq", ""))
                    if alloc1[index-1]:
                        out_fasta1.write(line)
                        flag = 1
                    else:
                        flag = 0
                else:
                    if flag:
                        out_fasta1.write(line)

    with open(out_prefix + '_reads_with_candidate_itrs_2.fasta', "w") as out_fasta2:
        flag = 0
        with open(read_fasta2, "r") as f2:
            for line in f2:
                if line[0] == ">":
                    index = int(line.split("_")[0].replace(">Seq", ""))
                    if alloc2[index-1]:
                        out_fasta2.write(line)
                        flag = 1
                    else:
                        flag = 0
                else:
                    if flag:
                        out_fasta2.write(line)


def get_contigs_with_itrs(sam_file, alloc_array, offset_array, read_pair, insert_size, read_length):

    # Get contigs and ITR cluster positions
    ir_positions = []
    ir_paired_reads = []
    itr_info = defaultdict(lambda: '')

    read_order = ['_f1', '_f2']
    if read_pair == 2:
        read_order = ['_f2', '_f1']

    with open(sam_file, "r") as sam:
        first_contig = ''

        for line in sam:
            current_contig = line.split("\t")[2]
            sam_flag = int(line.split('\t')[1])
            ir_paired_read = line.split("\t")[0]
            read_pair_id = ir_paired_read.split('_LCoord')[0][-3:]
            ir_position = 0
            index = int(ir_paired_read.split("_")[0].replace("Seq", ""))

            if read_pair_id == read_order[0]:
                if sam_flag in (99, 147, 83, 163, 81, 161, 97, 145, 65, 129, 113, 177):
                    read_position = int(line.split("\t")[3])
                    ir_offset = int(ir_paired_read.split('_RCoord')[0].split('LCoord_')[1])
                    ir_position = read_position + ir_offset - 1

            elif read_pair_id == read_order[1]:
                if sam_flag not in (99, 147, 83, 163, 81, 161, 97, 145, 65, 129, 113, 177):
                    read_position = int(line.split("\t")[3])
                    ir_position = read_position

            if first_contig == '' and ir_position:
                ir_positions.append(ir_position)
                ir_paired_reads.append(ir_paired_read)
                first_contig = current_contig

            elif first_contig == current_contig and ir_position:
                for ind, pos in enumerate(ir_positions[::-1]):
                    if abs(ir_position - pos) <= 2750:
                        if abs(ir_position - pos) >= 500:
                            index = int(ir_paired_reads[ind].split("_")[0].replace("Seq", ""))
                            alloc_array[index-1] = 1
                            if itr_info[current_contig] == '':
                                itr_info[current_contig] = defaultdict(lambda: '')
                                itr_info[current_contig][ir_paired_reads[ind]] = pos
                            else:
                                if itr_info[current_contig][ir_paired_reads[ind]] == '':
                                    itr_info[current_contig][ir_paired_reads[ind]] = pos
                            index = int(ir_paired_read.split("_")[0].replace("Seq", ""))
                            alloc_array[index-1] = 1
                            itr_info[current_contig][ir_paired_read] = ir_position
                    else:
                        if ind:
                            ir_positions_tmp = ir_positions[-ind:]
                            ir_paired_reads_tmp = ir_paired_reads[-ind:]
                            ir_positions = ir_positions_tmp
                            ir_paired_reads = ir_paired_reads_tmp
                        else:
                            ir_positions = []
                            ir_paired_reads = []
                        break

                ir_positions.append(ir_position)
                ir_paired_reads.append(ir_paired_read)

            elif ir_position:
                ir_positions = [ir_position]
                ir_paired_reads = [ir_paired_read]
                first_contig = current_contig

    return([itr_info, alloc_array])


def get_arguments():
    parser = argparse.ArgumentParser(description='Get reads that contain candidate ITRs and contigs that associate with these.')
    parser.add_argument('--contig_fasta', '-c', dest='contig_fasta', required=True,
                        help='Input contig FASTA file.', type = str)
    parser.add_argument('--sam_file1', '-s1', dest='sam_file1', required=True,
                        help='Input sam file from mapping of reads to contigs.', type = str)
    parser.add_argument('--sam_file2', '-s2', dest='sam_file2', required=True,
                        help='Input sam file from mapping of reads to contigs.', type = str)
    parser.add_argument('--fasta1', '-f1', dest='fasta_file1', required=True,
                        help='Input first paied FASTA file.', type = str)
    parser.add_argument('--fasta2', '-f2', dest='fasta_file2', required=True,
                        help='Input second paied FASTA file.', type = str)
    parser.add_argument('--insert_size', '-i', dest='insert_size', required=True,
                        help='Predicted insert size between paired reads.', type = float)
    parser.add_argument('--read_length', '-l', dest='read_length', required=True,
                        help='Length of reads.', type = float)
    parser.add_argument('--output_prefix', '-o', dest='output_prefix', required=True,
                        help='Prefix of output files.', type = str)
    return parser


def main(args):
    # Allocate for fasta indexing
    size = get_array_size(args.fasta_file1, args.fasta_file2)

    print("Allocating empty arrays...")
    alloc_arrays = allocate_seqs(size)

    print("Allocating empty arrays with IR size offset...")
    offset_arrays = get_ir_offset(args.fasta_file1, args.fasta_file1, size)

    # Get contigs, reads paired to reads with candidate ITRs and mapping positions
    print("Getting ITR positions...")
    contigs_read_and_positions1 = get_contigs_with_itrs(args.sam_file1, alloc_arrays[0], offset_arrays[0], 1, args.insert_size, args.read_length)
    contigs_read_and_positions2 = get_contigs_with_itrs(args.sam_file2, alloc_arrays[1], offset_arrays[1], 2, args.insert_size, args.read_length)

    # Write contig FASTA file
    print("Writing contig FASTA...")
    write_contig_file(args.contig_fasta, contigs_read_and_positions1[0], contigs_read_and_positions2[0], args.output_prefix)

    # Write FASTA files and get reads containing candidate ITRs
    print("Writing reads with ITRs FASTAs...")
    write_fasta_files(args.fasta_file1, args.fasta_file2, (contigs_read_and_positions1[1], contigs_read_and_positions2[1]), args.output_prefix)

    # Write tab file with contig and read names and mapping positions
    print("Writing tab file...")
    write_tab_file(contigs_read_and_positions1[0], contigs_read_and_positions2[0], args.output_prefix)


if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))
