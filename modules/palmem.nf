process palmem {

    input:
	tuple val(sample_id), path(fasta1), path(fasta2)

	output:
    tuple val(sample_id), path("${sample_id}_IR_1.fasta"), path("${sample_id}_IR_2.fasta"), path("${sample_id}_IR.tab"), emit: ir_ch

    script:
    min_itr_length = params.min_itr_length
    max_itr_length = params.max_itr_length
    kmer_length = params.kmer_length
	"""
	pal-mem -f1 ${fasta1} -f2 ${fasta2} -t ${task.cpus} -l ${min_itr_length} -m ${max_itr_length} -k ${kmer_length} -o ${sample_id}
	"""
}
