/*
 * Create IS catalog
 */
process createISCatalog {
    input:
    path clstr_file
    path itr_read_tab_files
    path is_annot_tab_files

    output:
    path "${output_prefix}_insertion_sequence_catalog.txt"

    script:
    output_prefix=params.batch_name

    """
    create_IS_catalog.py \
        --clstr_file ${clstr_file} \
        --itr_reads_tab ${itr_read_tab_files} \
        --is_annot_tab ${is_annot_tab_files} \
        --output_prefix ${output_prefix}
    """
}
