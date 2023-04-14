process getInsertionSequences {

    input:
    tuple val(sample_id), path(fasta), path(contig_info), path(tsv)

    output:
    path("${output_txt}"), emit: txt
    path("${output_fasta}"), emit: fasta

    script:
    contigs_out="${sample_id}_filtered_transposase.fasta"
    info_out="${sample_id}_filtered_transposase.tsv"
    output_txt="${sample_id}_insertion_sequences_info.txt"
    output_fasta="${sample_id}_insertion_sequences.fasta"

    """
    get_contigs_with_transposase.py -f ${fasta} -t ${tsv} -o ${sample_id}

    get_insertion_sequences.py \
        --contig_fasta ${contigs_out} \
        --contig_info ${contig_info} \
        --interproscan_info ${info_out} \
        --output_prefix ${sample_id}
    """
}
