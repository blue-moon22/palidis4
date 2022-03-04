/*
 * Get info on read annotations and their sequences for ITR clusters
 */
process getReadAndFastaInfo {

    input:
    tuple val(sample_id), file(is_info_tab), file(ir_info_tab), file(ir_fasta)

    output:
    path("${output_prefix}_ITRs.fasta"), emit: itr_fasta_ch
    path("${sample_id}_itr_read_positions_clusters.txt"), emit: itr_clusters_ch

    script:


    """
    get_read_and_fasta_info.py \
        --ir_fasta ${ir_fasta} \
        --ir_cluster_tab ${ir_info_tab} \
        --is_annotations_tab ${is_info_tab}
        --output_prefix ${sample_id}
    """
}
