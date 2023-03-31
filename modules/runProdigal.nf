/*
 * Run Prodigal
 */
process runProdigal {

    input:
    tuple val(sample_id), file(fasta)

    output:
    tuple val(sample_id), path("${sample_id}.faa"), emit: prot, optional: true
    tuple val(sample_id), path(fasta), emit: nucl, optional: true

    script:
    """
    if [ -s ${fasta} ]
    then
        prodigal -p meta -i ${fasta} -a ${sample_id}.faa -o ${sample_id}.gbk
    fi
    """
}
