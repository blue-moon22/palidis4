process mapReads {

	input:
	tuple val(sample_id), path(fasta1), path(fasta2), path(db_path), val(pair)

	output:
	tuple val(sample_id), path("${prefix}.sam.mapped.sorted")

    script:
    prefix = "${sample_id}_itr_in_${pair}"
	"""
    tar -xf ${db_path}
	bowtie2 --very-sensitive-local -x ${sample_id}_contigs -1 ${fasta1} -2 ${fasta2} -S ${prefix}.sam -p ${task.cpus} -f
	samtools view -S -b ${prefix}.sam -@ ${task.cpus} > ${prefix}.bam
	rm ${prefix}.sam
	samtools view -b -F 4 ${prefix}.bam -@ ${task.cpus} > ${prefix}.bam.mapped
	rm ${prefix}.bam
	samtools sort ${prefix}.bam.mapped -o ${prefix}.bam.mapped.sorted -@ ${task.cpus}
	rm ${prefix}.bam.mapped
	samtools index ${prefix}.bam.mapped.sorted
	samtools view ${prefix}.bam.mapped.sorted > ${prefix}.sam.mapped.sorted
	rm ${prefix}.bam.mapped.sorted*
    rm ${sample_id}_contigs*
	"""
}
