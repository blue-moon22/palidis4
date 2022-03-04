/*
 * Combine IR sequences
 */
process createIRCatalog {
    input:
    path ir_fastas

    output:
    path "all_ITRs.fasta"

    script:
    """
    cat ${ir_fastas} > all_IRs.fasta
    """
}
