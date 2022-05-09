/*
 * Get IS Information
 */
process getISInfoWithCOBS {

    container 'bluemoon222/palidis_dependencies:0.0.17'

    input:
    tuple val(sample_id), path(is_tab_file), path(isfinder_blast_out), path(cobs_out), file(isfinder_info_csv)

    output:
    path("${output}")

    script:
    output="${sample_id}_insertion_sequences_info.txt"
    """
    # Get organism from COBS bio sample ids
    ids=\$(sed '1d' ${cobs_out} | cut -f2 | sed -e :a -e '\$!N; s/\n/ /; ta')
    ffq \$ids > ${sample_id}_ffq.json

    get_IS_info.py --blast_out ${isfinder_blast_out} \
        --tab_file ${is_tab_file} \
        --is_finder_info ${isfinder_info_csv} \
        --cobs_search_out ${cobs_out} \
        --ffq_json ${sample_id}_ffq.json \
        --output_prefix ${sample_id}
    """
}

process getISInfoWithoutCOBS {

    container 'bluemoon222/palidis_dependencies:0.0.16'

    input:
    tuple val(sample_id), path(is_tab_file), path(isfinder_blast_out), file(isfinder_info_csv)

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
