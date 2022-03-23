/*
 * Combine ITR sequences
 */
process createITRCatalog {
    input:
    path itr_fastas

    output:
    path "all_ITRs.fasta"

    script:
    """
    cat ${itr_fastas} > all_ITRs.fasta
    """
}
