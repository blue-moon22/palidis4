/*
 * Create IS catalog
 */
process createISCatalog {
    input:
    path clstr_file
    path txt_files
    path tab_files

    output:
    path "${output_prefix}_insertion_sequence_annotations_catalog.tab"

    script:
    output_prefix=params.batch_name

    """
    process_clstr_file.py -c ${clstr_file} -o all.clstr
    assign_ITR_clusters.R -c all.clstr.tab -t _reads_itr_clusters.txt -a _insertion_sequence_annotations.tab -o ${output_prefix}
    """
}
