/*
 * Convert FASTQ.GZ to FASTA, replace spaces in names with underscores, add f1/f2 and Seq number
 */
process convertToFasta {

	input:
    tuple val(sample_id), path(reads_folder)

	output:
	tuple val(sample_id), path("${sample_id}.fasta")

    script:
    fastq = "${sample_id}.fastq"

	"""
    for file in ${reads_folder}/*fastq.gz
    do
        gunzip -c \$file >> ${fastq}
    done

    convert_fastq_to_fasta.py -f ${fastq} -o ${sample_id}.fasta

    rm ${fastq}
    """
}
