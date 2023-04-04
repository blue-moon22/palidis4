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
    num_lines=\$(echo \$(( \$(cat ${fasta} | wc -l) / 2 )))
    head -n \$num_lines ${fasta} > ${fasta}1
    tail -n \$num_lines ${fasta} > ${fasta}2

	pal-mem -f1 ${fasta}1 -f2 ${fasta}2 -t ${task.cpus} -l ${min_itr_length} -m ${max_itr_length} -k ${kmer_length} -o ${sample_id}
    cat ${sample_id}_IR_1.fasta ${sample_id}_IR_2.fasta > ${sample_id}_IR.fasta

    rm ${fasta}1 ${fasta}2 ${sample_id}_IR_1.fasta ${sample_id}_IR_2.fasta
	"""
}
