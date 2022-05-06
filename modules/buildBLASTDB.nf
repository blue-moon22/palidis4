process buildBLASTDB {

    input:
    file(fasta)

    output:
    path("blast_db")

    script:
    """
    mkdir blast_db
    mv ${fasta} blast_db/
    makeblastdb -in blast_db/${fasta} -out blast_db/${fasta} -parse_seqids -dbtype nucl
    """
}
