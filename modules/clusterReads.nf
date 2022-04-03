/*
 * Cluster reads
 */
process clusterReads {

    input:
    tuple val(sample_id), file(read_file)

    output:
    tuple val(sample_id), path("${output_prefix}.fasta.clstr"), emit: cluster_ch
    tuple val(sample_id), path("${sample_id}_irs.fasta"), emit: clipped_read_ch

    script:
    G=params.cd_hit_G
    aL=params.cd_hit_aL
    aS=params.cd_hit_aS
    A=params.min_itr_length
    c=params.cd_hit_c
    output_prefix="${sample_id}_nonred_G${G}_aL${aL}_aS${aS}_A${A}"

    """
    clip_reads.py --read_fasta ${read_file} --output_prefix ${sample_id}
    cd-hit-est -i ${sample_id}_irs.fasta -o ${output_prefix}.fasta -c ${c} -G ${G} -aL ${aL} -aS ${aS} -A ${A} -M 64000 -T ${task.cpus} -d 0
    """
}
