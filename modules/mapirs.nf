process mapIRs {

	input:
	tuple val(sample_id), path(fasta), path(db_path)

	output:
	tuple val(sample_id), path("${sample_id}.sam.mapped.sorted")

    script:
	"""
    tar -xf ${db_path}

	bowtie2 --very-sensitive-local -x ${sample_id}_contigs -U ${fasta} -S ${sample_id}.sam -p ${task.cpus} -f
	rm ${sample_id}_contigs*

	samtools view -S -b ${sample_id}.sam -@ ${task.cpus} > ${sample_id}.bam
	rm ${sample_id}.sam

	samtools view -b -F 4 ${sample_id}.bam -@ ${task.cpus} > ${sample_id}.bam.mapped
	rm ${sample_id}.bam

	samtools sort ${sample_id}.bam.mapped -o ${sample_id}.bam.mapped.sorted -@ ${task.cpus}
	rm ${sample_id}.bam.mapped

	samtools index ${sample_id}.bam.mapped.sorted
	samtools view ${sample_id}.bam.mapped.sorted > ${sample_id}.sam.mapped.sorted
	rm ${sample_id}.bam.mapped.sorted*
	"""
}
