#!/usr/bin/env python3
import sys, getopt
import pickle
from collections import defaultdict

def get_contigs_with_itrs(sam_file, ir_index_file):

    # Read ITR index file
    with open(ir_index_file, "rb") as input:
        itr_index = pickle.load(input)

    # Get contigs and ITR cluster positions
    itr_cluster_positions = defaultdict(lambda: '')
    with open(sam_file, "r") as sam:
        contigs_itr_clusters = defaultdict(lambda: '')
        first_contig = ''
        for line in sam:
            current_contig = line.split("\t")[2]
            if current_contig == first_contig:
                itr_cluster = itr_index[int(line.split("Seq")[1].split("_")[0]) - 1][line.split(".1_")[0].split(".2_")[0].split("_")[1].split("/1")[0].split("/2")[0]]
                if itr_cluster_positions[itr_cluster] == '':
                    itr_cluster_positions[itr_cluster] = [int(line.split("\t")[3])]
                else:
                    tmp = itr_cluster_positions[itr_cluster]
                    tmp.append(int(line.split("\t")[3]))
                    itr_cluster_positions[itr_cluster] = tmp
            else:
                for itr_cluster in itr_cluster_positions:
                    posx = itr_cluster_positions[itr_cluster]
                    if len(posx) > 1:
                        # Only accept ITR clusters that are positioned appropriately sequentially
                        pos_diff = [posx[n]-posx[n-1] for n in range(1,len(posx))]
                        itr_positions = [(posx[:-1][i], posx[:-1][i] + n) for i, n in enumerate(pos_diff) if n >= 500 and n <= 2750]
                        if itr_positions:
                            if contigs_itr_clusters[current_contig] == '':
                                contigs_itr_clusters[current_contig] = [str(itr_cluster) + '\t' + ';'.join(str(itr_positions))]
                            else:
                                tmp = contigs_itr_clusters[current_contig]
                                contigs_itr_clusters[current_contig] = tmp.append(str(itr_cluster) + '\t' + ';'.join(str(itr_positions)))
                # Reset
                itr_clusters_positions = defaultdict(lambda: '')
                first_contig = current_contig
    return(contigs_itr_clusters)


def write_files(contigs_dict, contigs_fasta_file, out_prefix):
    # Write fasta file
    with open(out_prefix + '_contigs_with_itrs.fa', "w") as out_fasta:
        flag = 0
        with open(contigs_fasta_file, "r") as f:
            for line in f:
                if line[0] == ">":
                    if line.split("\n")[0].replace(">", "") in contigs_dict.keys():
                        out_fasta.write(line)
                        flag = 1
                    else:
                        flag = 0
                else:
                    if flag:
                        out_fasta.write(line)
    # Write tab file
    with open(out_prefix + '_contigs_itr_info.tab', "w") as out_tab:
        out_tab.write('contig\titr_cluster\titr_positions\n')
        for contig in contigs_dict:
            for info in contigs_dict[contig]:
                out_tab.write(contig + '\t' + info[0] + '\n')


def main(argv):
    try:
        opts, args = getopt.getopt(argv, "h")
    except getopt.GetoptError:
        print("get_contigs_with_sites.py <contig_fasta> <contig_mapping_sam_file> <ir_index_file> <out_prefix>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print("get_contigs_with_sites.py <contig_fasta> <contig_mapping_sam_file> <ir_index_file> <out_prefix>")
            sys.exit()
    fasta = sys.argv[1]
    sam_file = sys.argv[2]
    itr_index_file = sys.argv[3]
    out_prefix = sys.argv[4]

    contigs_with_itr_positions = get_contigs_with_itrs(sam_file, itr_index_file)
    write_files(contigs_with_itr_positions, fasta, out_prefix)


if __name__ == "__main__":
    main(sys.argv[1:])
