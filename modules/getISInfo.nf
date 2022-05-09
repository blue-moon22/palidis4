/*
 * Get IS Information
 */
process getISInfoWithCOBS {

    input:
    tuple val(sample_id), path(is_tab_file), path(isfinder_blast_out), path(cobs_out), file(isfinder_info_csv)

    output:
    path("${output}")

    script:
    output="${sample_id}_insertion_sequences_info.txt"
    """
    set +e

    grep SAMN ${cobs_out} | cut -f2 | sort | uniq > SAMN_ids.txt
    num_ids=\$(cat SAMN_ids.txt | wc -l)

    for ((i=1;i<=\${num_ids};i++))
    do
        biosample_id=\$(sed -n "\${i}p" SAMN_ids.txt)
        ffq \$biosample_id > \${biosample_id}.json
    done

    get_IS_info.py --blast_out ${isfinder_blast_out} \
        --tab_file ${is_tab_file} \
        --is_finder_info ${isfinder_info_csv} \
        --cobs_search_out ${cobs_out} \
        --ffq_json \$(pwd) \
        --output_prefix ${sample_id}
    """
}

process getISInfoWithoutCOBS {

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
