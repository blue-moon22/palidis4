process searchISfinder {

    input:
    tuple val(sample_id), file(fasta), path(blast_db)

    output:
    tuple val(sample_id), path("${output}")

    when:
    fasta.size() > 0

    script:
    e_value=params.e_value
    isfinder_seq_fna=params.isfinder_seq_fna
    output="${sample_id}_blast.out"

    """
    blastn -query ${fasta} -subject ${blast_db}/${isfinder_seq_fna} -outfmt 6 -evalue ${e_value} > ${output}
    """
}
