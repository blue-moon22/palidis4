process buildBLASTDB {

    input:
    file(fasta)

    output:
    path("blast_db")

    script:
    """
    makeblastdb -in ${fasta} -out blast_db/${fasta} -parse_seqids -dbtype nucl
    """
}
