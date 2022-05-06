process searchISfinder {

    input:
    tuple val(sample_id), file(fasta), file(blast_db)

    output:
    tuple val(sample_id), path("${output}")

    script:
    e_value=params.e_value
    output="${sample_id}_blast.out"

    """
    blastn -query ${fasta} -subject ${blast_db} -outfmt 6 -evalue ${e_value} > ${output}
    """
}
