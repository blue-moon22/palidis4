/*
 * Get contigs that have transposase
 */
process contigCandidates {
    input:
    tuple val(sample_id), path(fasta), path(tsv)

    output:
    tuple val(sample_id), path(contigs_out), emit: ref
    tuple val(sample_id), path(contigs_out), path(info_out), emit: fasta_info

    script:
    contigs_out="${sample_id}_filtered_transposase.fasta"
    info_out="${sample_id}_filtered_transposase.tsv"

    """
    get_contigs_with_transposase.py -f ${fasta} -t ${tsv} -o ${sample_id}
    """
}
