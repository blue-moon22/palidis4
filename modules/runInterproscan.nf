/*
 * Run Interproscan
 */
process runInterproscan {

    input:
    tuple val(sample_id), path(faa), path(db)

    output:
    tuple val(sample_id), path(faa), path("${faa}.tsv"), optional: true

    script:
    lsf=params.lsf
    """
    if [ -s ${faa} ]
    then
        # Remove * from protein prediction
        sed -i 's/*//' ${faa}
        if ${lsf}
        then
            ./${db}/interproscan.sh -mode cluster -clusterrunid ${sample_id}_interproscan -i ${faa} -f tsv -dp -cpu ${task.cpus}
        else
            ./${db}/interproscan.sh -i ${faa} -f tsv -dp -cpu ${task.cpus}
        fi
    fi
    """
}