/*
 * Get IS Information
 */
process getISInfoWithCOBS {

    input:
    tuple val(sample_id), path(is_tab_file), path(isfinder_blast_out), path(cobs_out)
    file(isfinder_info_csv)

    output:
    path("${output}")

    script:
    output = "${sample_id}_insertion_sequences_info.txt"
    """
    get_IS_info.py --blast_out ${isfinder_blast_out} \
        --tab_file ${is_tab_file} \
        --is_finder_info ${isfinder_info_csv} \
        --cobs_search_out ${cobs_out} \
        --output_prefix ${sample_id}
    """
}

process getISInfoWithoutCOBS {

    input:
    tuple val(sample_id), path(is_tab_file), path(isfinder_blast_out)
    file(isfinder_info_csv)

    output:
    path("${output}")

    script:
    output = "${sample_id}_insertion_sequences_info.txt"
    """
    get_IS_info.py --blast_out ${isfinder_blast_out} \
        --tab_file ${is_tab_file} \
        --is_finder_info ${isfinder_info_csv} \
        --output_prefix ${sample_id}
    """
}
