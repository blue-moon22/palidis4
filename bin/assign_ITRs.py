#!/usr/bin/env python3
import argparse, sys, os
import subprocess
import multiprocessing as mp

from collections import defaultdict

def create_assembly_bins(assembly_file):

    assembly_bins_dict = {}
    with open(assembly_file, "r") as fa:
        for line in fa:
            if line[0] == '>':
                header = line[1:].replace('\n', '')
            else:
                assembly_bins_dict[header] = line.replace('\n', '')

    return assembly_bins_dict


def create_cluster_dictionary(itr_clusters):

    seq_num = 0
    with open(itr_clusters, "r") as fa:
        for line in fa:
            if line[0] != ">":
                curr_seq_num = int(line.split(">Seq")[1].split("_")[0])
                if curr_seq_num > seq_num:
                    seq_num = curr_seq_num
    alloc1 = ['']*seq_num
    alloc2 = ['']*seq_num

    with open(itr_clusters, "r") as cl:
        for line in cl:
            if line[0] == ">":
                cluster = line.split(" ")[1].split("\n")[0]
            else:
                seq_no = int(line.split(">Seq")[1].split("_")[0]) - 1
                if '_f1' in line:
                    alloc1[seq_no] = cluster
                elif '_f2' in line:
                    alloc2[seq_no] = cluster
    return((alloc1, alloc2))


def bin_positions(cl_dict, tab_info_file, assembly_bins_dict, output_prefix):

    clusters_positions = {}

    sample_id = output_prefix.split('/')[len(output_prefix.split('/'))-1]
    with open(output_prefix + '_reads_itr_clusters.txt', 'w') as out:
        out.write('sample_id\tcontig\tread\tir_start\tir_end\titr_cluster\tir\n')
        with open(tab_info_file, "r") as tab:
            next(tab)
            for line in tab:
                seq = line.split('\t')[2]
                start_pos = int(line.split('\t')[4])
                end_pos = int(line.split('\t')[5])
                seq_no = int(seq.replace('Seq', '').split('_')[0]) - 1
                current_contig = line.split('\t')[1]
                ir = line.split('\t')[3]
                if '_f1' in seq:
                    if ir:
                        cluster = cl_dict[0][seq_no]
                    else:
                        cluster = cl_dict[1][seq_no]
                elif '_f2' in seq:
                    if ir:
                        cluster = cl_dict[1][seq_no]
                    else:
                        cluster = cl_dict[0][seq_no]

                out.write(sample_id + '\t' + current_contig + '\t' + seq + '\t' + str(start_pos) + '\t' + str(end_pos) + '\t' + cluster + '\t' + str(ir) + '\n')

                if ir:

                    if current_contig not in clusters_positions.keys():
                        clusters_positions[current_contig] = {}

                    if cluster in clusters_positions[current_contig].keys():
                        encode_contig = list(clusters_positions[current_contig][cluster])
                    else:
                        encode_contig = ['0']*len(assembly_bins_dict[current_contig])

                    # Update positions
                    contig_len = len(assembly_bins_dict[current_contig])
                    if start_pos < 1:
                        if end_pos > contig_len:
                            encode_contig = ['1']*contig_len
                        else:
                            encode_contig[:end_pos] = ['1']*end_pos
                    else:
                        if end_pos > contig_len:
                            encode_contig[start_pos-1:contig_len] = ['1']*(contig_len-start_pos+1)
                        else:
                            encode_contig[start_pos-1:end_pos] = ['1']*(end_pos-start_pos+1)

                    clusters_positions[current_contig][cluster] = ''.join(encode_contig)

    return clusters_positions


def count_bins(sbin):
    count = 1
    bin = list(sbin)
    bin.append('NA')
    val = bin[0]
    count_bin_out = []
    for new_val in bin[1:]:
        if new_val == val:
            count += 1
        else:
            count_bin_out.append((val, count))
            val = new_val
            count = 1

    return count_bin_out


def get_itrs_from_count_bins(count_bins, MIN_IS_LEN, MAX_IS_LEN, MIN_ITR_LEN, MAX_ITR_LEN):
    is_gaps = []

    # Remove Ns from count bin and add 0s. Record position and Ns
    count_bins_mod = []
    sum_0s = 0
    sum_Ns = 0
    N_bin = ()
    N_pos = 0
    for ind, count_bin in enumerate(count_bins):
        if count_bin[0] != 'N':
            if ind > 0 and ind < len(count_bins)-1:
                if count_bins[ind+1][0] == 'N':
                    N_pos = ind
                    sum_0s += count_bin[1]
                elif count_bins[ind-1][0] == 'N':
                    sum_0s += count_bin[1]
                    N_bin = (N_pos, sum_Ns)
                    count_bins_mod.append(('0',sum_0s))

                else:
                    count_bins_mod.append(count_bin)

            else:
                count_bins_mod.append(count_bin)

        else:
            sum_Ns += count_bin[1]

    count_bins = count_bins_mod

    for ind, count_bin in enumerate(count_bins):
        if count_bin[0] == '0':
            if count_bin[1] >= MIN_IS_LEN and count_bin[1] <= MAX_IS_LEN:
                if ind != 0 and ind != len(count_bins) - 1:
                    is_gaps.append(ind)

    itr_positions = []
    for ind in is_gaps:
        if count_bins[ind-1][1] >= MIN_ITR_LEN and count_bins[ind-1][1] <= MAX_ITR_LEN and count_bins[ind+1][1] >= MIN_ITR_LEN and count_bins[ind+1][1] <= MAX_ITR_LEN:
            itr1_start = sum(e[1] for e in count_bins[:ind-1]) + 1
            itr1_end = sum(e[1] for e in count_bins[:ind-1]) + count_bins[ind-1][1]
            itr2_start = sum(e[1] for e in count_bins[:ind+1]) + 1
            itr2_end = sum(e[1] for e in count_bins[:ind+1]) + count_bins[ind+1][1]

            if N_bin:
                if ind >= N_bin[0]:
                    itr2_start += N_bin[1]
                    itr2_end += N_bin[1]
                    if ind > N_bin[0]:
                        itr1_start += N_bin[1]
                        itr1_end += N_bin[1]

            itr_positions.append((itr1_start, itr1_end, itr2_start, itr2_end))

    return itr_positions


def get_itr_sequences(contig_seq, positions):

    itr1 = contig_seq[positions[0]-1:positions[1]]
    itr2 = contig_seq[positions[2]-1:positions[3]]

    return [itr1, itr2]


def check_blast_out(out_blast_file, MIN_ITR_LEN):
    flag = 0
    with open(out_blast_file) as tmp:
        for line in tmp: # Check orientation (most complete) match
            if "Identities =" in line:
                local_len = int(line.split(" ")[3].split("/")[1])
                if local_len >= MIN_ITR_LEN:
                    flag = 1
                else:
                    flag = 0
            elif "Strand=" in line:
                if flag:
                    if "Plus/Minus" in line or "Minus/Plus" in line:
                        return True
                    else:
                        return False


def are_reverse_cmp(itr_sequences, MIN_ITR_LEN):
    result = 0
    # Write temporary FASTAs
    itr1_fasta = 'itr1_tmp.fasta'
    with open(itr1_fasta, 'w') as out:
        out.write('>itr1\n' + itr_sequences[0] + '\n')

    itr2_fasta = 'itr2_tmp.fasta'
    with open(itr2_fasta, 'w') as out:
        out.write('>itr2\n' + itr_sequences[1] + '\n')

    # Run blastn
    out_blast_file = "blastn_{}_{}.out".format(itr_sequences[0], itr_sequences[1])
    run_blastn_cmd = "blastn -query {} -subject {} -task blastn -word_size 4 > {}".format(itr1_fasta, itr2_fasta, out_blast_file)
    os.system(run_blastn_cmd)

    return check_blast_out(out_blast_file, MIN_ITR_LEN)


def remove_clusters_positions(clusters, itr_positions):
    clusters_new = {}
    for cluster, bin in clusters.items():
        for pos in itr_positions:
            bin = list(bin)
            bin[pos[0]-1:pos[3]] = ['N']*(pos[3] - (pos[0]-1))
            bin = ''.join(bin)
        clusters_new[cluster] = bin

    return clusters_new


def annotate_itrs(clusters, contig_seq, MIN_IS_LEN, MAX_IS_LEN, MIN_ITR_LEN, MAX_ITR_LEN, output_info):
    itr_clusters = {}
    for cluster, bin in clusters.items():
        count_bins_out = count_bins(bin)
        itr_positions = get_itrs_from_count_bins(count_bins_out, MIN_IS_LEN, MAX_IS_LEN, MIN_ITR_LEN, MAX_ITR_LEN)
        for pos in itr_positions:
            itr_sequences = get_itr_sequences(contig_seq, pos)
            print(itr_sequences)
            if are_reverse_cmp(itr_sequences, MIN_ITR_LEN):
                if pos in itr_clusters:
                    itr_pos = itr_clusters[cluster]
                    itr_pos.append(pos)
                    itr_clusters[cluster] = itr_pos
                else:
                    itr_clusters[cluster] = [pos]
                output_info.append(str(pos[0]) + '\t' + str(pos[1]) + '\t' + str(pos[2]) + '\t' + str(pos[3]) + '\t' + cluster + '\n')

    # Rerun function with updated clusters_positions to find nested ITRs
    for cluster, positions in itr_clusters.items():
        clusters_new = remove_clusters_positions(clusters, positions)
        annotate_itrs(clusters_new, contig_seq, MIN_IS_LEN, MAX_IS_LEN, MIN_ITR_LEN, MAX_ITR_LEN, output_info)

    return output_info


def write_itr_annotations(clusters_positions, assemblies_dict, MIN_IS_LEN, MAX_IS_LEN, MIN_ITR_LEN, MAX_ITR_LEN, output_prefix):

    clusters_itrs = {}
    sample_id = output_prefix.split('/')[len(output_prefix.split('/'))-1]

    pool = mp.Pool(mp.cpu_count())
    output = [pool.apply(annotate_itrs, args=(clusters, assemblies_dict[contig], MIN_IS_LEN, MAX_IS_LEN, MIN_ITR_LEN, MAX_ITR_LEN, [])) for contig, clusters in clusters_positions.items()]
    pool.close()

    with open(output_prefix + '_insertion_sequence_annotations.tab', 'w') as out:
        out.write('sample_id\tcontig\titr1_start_position\titr1_end_position\titr2_start_position\titr2_end_position\titr_cluster\n')
        for ind, positions_cluster in enumerate(output):
            for pos in positions_cluster:
                out.write(sample_id + '\t' + list(clusters_positions)[ind] + '\t' + pos)


def get_arguments():
    parser = argparse.ArgumentParser(description='Get the ITR clusters and information on contigs and reads with ITRs.')
    parser.add_argument('--cdhit_cluster_file', '-c', dest='cluster_file', required=True,
                        help='Input .clstr file from CD-HIT.', type = str)
    parser.add_argument('--info_tab_file', '-t', dest='tab_file', required=True,
                        help='Input "contigs_reads_ir_position_info.tab" file.', type = str)
    parser.add_argument('--assemblies_fasta_file', '-a', dest='assemblies_file', required=True,
                        help='Input assemblies (one-line) FASTA file.', type = str)
    parser.add_argument('--min_is_len', '-min_is', dest='min_is_len', required=True,
                        help='Minimum length of insertion sequence', type = int, default = 500)
    parser.add_argument('--max_is_len', '-max_is', dest='max_is_len', required=True,
                        help='Maximum length of insertion sequence', type = int, default = 3000)
    parser.add_argument('--min_itr_len', '-min_itr', dest='min_itr_len', required=True,
                        help='Minimum length of insertion sequence', type = int, default = 25)
    parser.add_argument('--max_itr_len', '-max_itr', dest='max_itr_len', required=True,
                        help='Maximum length of insertion sequence', type = int, default = 50)
    parser.add_argument('--output_prefix', '-o', dest='output_prefix', required=True,
                    help='Prefix of output files.', type = str)
    return parser


def main(args):

    # Set min and max IS length
    MAX_IS_LEN = args.max_is_len
    MIN_IS_LEN = args.min_is_len
    MAX_ITR_LEN = args.max_itr_len
    MIN_ITR_LEN = args.min_itr_len

    # Create a dictionary of assemblies with binary values of assembly length
    assembly_bins_dict = create_assembly_bins(args.assemblies_file)

    # Create dictionary of IR clusters
    cl_dict = create_cluster_dictionary(args.cluster_file)

    # Populate binary values with 1s to show IR coverage and write tab file with read positions of IRs and their clusters
    clusters_positions = bin_positions(cl_dict, args.tab_file, assembly_bins_dict, args.output_prefix)

    # Write putative insertion sequences
    write_itr_annotations(clusters_positions, assembly_bins_dict, MIN_IS_LEN, MAX_IS_LEN, MIN_ITR_LEN, MAX_ITR_LEN, args.output_prefix)


if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))
