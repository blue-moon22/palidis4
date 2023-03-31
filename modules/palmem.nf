process palmem {

    input:
	tuple val(sample_id), path(fasta)

	output:
    tuple val(sample_id), path("${sample_id}_IR.fasta"), path("${sample_id}_IR.tab")

    script:
    min_itr_length = params.min_itr_length
    max_itr_length = params.max_itr_length
    kmer_length = params.kmer_length
	"""
	pal-mem -fu ${fasta} -t ${task.cpus} -l ${min_itr_length} -m ${max_itr_length} -k ${kmer_length} -o ${sample_id}
	"""
}
