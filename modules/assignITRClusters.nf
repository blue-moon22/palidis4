/*
 * Run CD-HIT on IRs
 */
process assignITRClusters {
    input:
	path all_itrs

	output:
    path "all.fasta.clstr"

    script:
    G=params.cd_hit_G
    aL=params.cd_hit_aL
    aS=params.cd_hit_aS
    A=params.min_itr_length
    output_prefix="all"

    """
    cd-hit-est -i ${all_itrs} -o ${output_prefix}.fasta -c 1.0 -G ${G} -aL ${aL} -aS ${aS} -A ${A} -M 64000 -T ${task.cpus} -d 0
    """
}
