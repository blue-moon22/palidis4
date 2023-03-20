/*
 * Run Interproscan
 */
process runInterproscan {

    input:
    tuple val(sample_id), path(faa), path(fasta), path(db)

    output:
    tuple val(sample_id), path(faa), path("${faa}.tsv"), emit: faa, optional: true
    tuple val(sample_id), path(fasta), path("${faa}.tsv"), emit: fasta, optional: true

    script:
    """
    # Remove * from protein prediction
    sed -i 's/*//' ${faa}

    ./${db}/interproscan.sh -i ${faa} -f tsv -dp -cpu ${task.cpus}
    """
}
