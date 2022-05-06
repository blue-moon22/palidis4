process buildBLASTDB {

    input:
    file(fasta)

    output:
    path("${fasta}.*")

    script:
    """
    makeblastdb -in ${fasta} -out ${fasta} -parse_seqids -dbtype nucl
    """
}
