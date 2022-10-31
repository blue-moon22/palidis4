/*
 * Convert FASTQ.GZ to FASTA, replace spaces in names with underscores, add f1/f2 and Seq number
 */
process convertToFasta {

	input:
    tuple val(sample_id), path(read1), path(read2)

	output:
	tuple val(sample_id), path("${sample_id}_1.fasta"), path("${sample_id}_2.fasta")

    script:
    fastq1 = "${sample_id}_1.fastq"
    fastq2 = "${sample_id}_2.fastq"
    chunk_size=params.chunk_size

	"""
    for file in ${read1}
    do
        gunzip -c \$file >> ${fastq1}
    done

    for file in ${read2}
    do
        gunzip -c \$file >> ${fastq2}
    done

    sed -n '1~4s/^@/>/p;2~4p' ${fastq1} > ${sample_id}_1.fa
    sed -n '1~4s/^@/>/p;2~4p' ${fastq2} > ${sample_id}_2.fa
    rm ${fastq1} ${fastq2}

    modify_headers.py -f ${sample_id}_1.fa -r 1
    modify_headers.py -f ${sample_id}_1.fa -r 2
    rm ${sample_id}_*.fa

    """
}
