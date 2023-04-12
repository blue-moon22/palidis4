/*
 * Concatenate files
 */
process mergeTSV {

    input:
    tuple val(sample_id), path(tsv)

    output:
    tuple val(sample_id), path(merged_file)

    script:
    merged_file="${sample_id}_merged.tsv"
    """
    cat ${tsv} >> ${merged_file}
    """
}
