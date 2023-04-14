/*
 * Filter contigs less than minimum IS length
 */
process filterContigs {
    input:
    tuple val(sample_id), path(contigs_path)

    output:
    tuple val(sample_id), path(contigs_out), emit: db_ch
    tuple val(sample_id), path(contigs_out), emit: fasta_ch

    script:
    min_is_len = params.min_is_len
    contigs_out = "assemblies_filtered.fasta"

    """
    remove_contigs_smaller_than.py -l ${min_is_len} -i ${contigs_path} -o ${contigs_out}
    """
}
