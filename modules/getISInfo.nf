/*
 * Get IS Information
 */
process getISInfoWithCOBS {

    input:
    tuple val(sample_id), path(is_tab_file), path(cobs_out), path(faa_file), path(tsv_file)

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

    get_IS_info.py \
        --tab_file ${is_tab_file} \
        --cobs_search_out ${cobs_out} \
        --ffq_json \$(pwd) \
        --aa_fasta ${faa_file} \
        --interproscan_out ${tsv_file} \
        --output_prefix ${sample_id}
    """
}

process getISInfoWithoutCOBS {

    input:
    tuple val(sample_id), path(is_tab_file), path(faa_file), path(tsv_file)

    output:
    path("${output}")

    script:
    output = "${sample_id}_insertion_sequences_info.txt"
    """
    get_IS_info.py \
        --tab_file ${is_tab_file} \
        --aa_fasta ${faa_file} \
        --interproscan_out ${tsv_file} \
        --output_prefix ${sample_id}
    """
}
