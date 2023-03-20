/*
 * Clip reads and cluster inverted repeats
 */
process clusterIRs {

    input:
    tuple val(sample_id), file(ir_1), file(ir_2)

    output:
    tuple val(sample_id), path(output_fasta), emit: clipped_ir_ch
    tuple val(sample_id), path("${output_prefix}.fasta.clstr"), emit: cluster_ch

    script:
    G=params.cd_hit_G
    aL=params.cd_hit_aL
    aS=params.cd_hit_aS
    A=params.min_itr_length
    c=params.cd_hit_c
    output_fasta=""${sample_id}_clipped_irs.fasta"
    output_prefix="${sample_id}_nonred_G${G}_aL${aL}_aS${aS}_A${A}"

    """
    cat ${ir_1} ${ir_2} > ${sample_id}_IR.fasta
    clip_reads.py --read_fasta ${sample_id}_IR.fasta --output_prefix ${sample_id}_clipped
    cd-hit-est -i ${sample_id}_clipped_irs.fasta -o ${output_prefix}.fasta -c ${c} -G ${G} -aL ${aL} -aS ${aS} -A ${A} -M 64000 -T ${task.cpus} -d 0
    """
}
