/*
 * Run Prodigal
 */
process runProdigal {

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}.faa")

    script:
    """
    prodigal -p meta -i ${fasta} -a ${sample_id}.faa -o ${sample_id}.gbk
    """
}
