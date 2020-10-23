/*
 * Nextflow pipeline for identifying insertion sequences from metagenomic data
 *
 * Author:
 * Victoria Carr <victoriacarr018@gmail.com
 *
 */


/*
 * Define parameters
 */
params.reads=""

params.outpath="$baseDir"
contig_result_path = file(params.outpath)

params.hmm = "$baseDir/db/Pfam-A_transposase.hmm"
hmm_name = file(params.hmm).name
hmm_dir = file(params.hmm).parent

/*
 * Create the 'read_pairs_ch' channel that emits tuples containing three elements:
 * the pair ID, the first read-pair file and the second read-pair file
 */
Channel
	.fromFilePairs( params.reads )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set { read_pairs_ch }

/*
 * Preprocess: Decompress fastq file, convert files, replace spaces in names with underscores, add f1/f2 and Seq number
 */
process convertToFasta {
	input:
  tuple val(pair_id), file(reads) from read_pairs_ch

	output:
	tuple val(pair_id), file("${pair_id}_1.fastq"), file("${pair_id}_2.fastq"), file("${pair_id}_1.fasta"), file("${pair_id}_2.fasta") into fasta_ch

	"""
	gunzip -c ${reads[0]} > ${pair_id}_1.fastq
	gunzip -c ${reads[1]} > ${pair_id}_2.fastq
  convertfastq2fasta.sh ${pair_id}_1.fastq ${pair_id}_2.fastq ${pair_id}
  """
}

/*
 * Step 1: Assemble reads
 */
process assembly {
	label 'spades'

	input:
	tuple val(pair_id), file(fastq1), file(fastq2), file(fasta1), file(fasta2) from fasta_ch

	output:
	tuple val(pair_id), file(fasta1), file(fasta2), file("${pair_id}_contigs.fasta") into assembly_ch

  """
  spades.py -1 ${fastq1} -2 ${fastq2} -t $task.cpus -k 21,33,55 --only-assembler --meta -o ${pair_id}
  mv ${pair_id}/contigs.fasta ${pair_id}_contigs.fasta
  """
}

/*
 * Step 2: Run MMSeq2 lincluster
 */
process mmseqs2 {
	label 'mmseqs2'

	input:
	tuple val(pair_id), file(fasta1), file(fasta2), file(contigs) from assembly_ch

	output:
	tuple val(pair_id), file(fasta1), file(fasta2), file("${pair_id}_DB_rep.fasta"), file("${pair_id}_DB_clu_seq.fasta"), file(contigs) into mmseq_ch

	"""
  echo 'Building database...'
	mmseqs createdb ${fasta1} ${fasta2} ${pair_id}_DB

	echo 'Running linclust...'
	mmseqs linclust ${pair_id}_DB ${pair_id}_DB_clu tmp --cov-mode 2 -c 0.5 --threads ${task.cpus}

	echo 'Convert sequences to fasta...'
	mmseqs createsubdb ${pair_id}_DB_clu ${pair_id}_DB ${pair_id}_DB_rep
	mmseqs convert2fasta ${pair_id}_DB_rep ${pair_id}_DB_rep.fasta

	echo 'Generate TSV output...'
	mmseqs createtsv ${pair_id}_DB ${pair_id}_DB ${pair_id}_DB_clu ${pair_id}_DB_clu.tsv

	echo 'Get sequence information...'
	mmseqs createseqfiledb ${pair_id}_DB ${pair_id}_DB_clu ${pair_id}_DB_clu_seq
	mmseqs result2flat ${pair_id}_DB ${pair_id}_DB ${pair_id}_DB_clu_seq ${pair_id}_DB_clu_seq.fasta
	"""
}

/*
 * Step 3: Run pal-MEM
 */
process palmem {
	label 'pal_mem'

	input:
	tuple val(pair_id), file(fasta1), file(fasta2),	file(rep_fasta), file(clu_fasta), file(contigs) from mmseq_ch

	output:
	tuple val(pair_id), file(fasta1), file(fasta2), file("${pair_id}_palmem_ITR.fasta"), file(contigs) into palmem_ch

	"""
	pal-mem -fu ${rep_fasta} -t ${task.cpus} -l ${params.mem_length} -k ${params.kmer_length} -o ${pair_id}_palmem
	"""
}

/*
 * Step 4: Get discordant reads
 */
process getDiscordantReads {
  label 'discord_reads'

	input:
	tuple val(pair_id), file(fasta1), file(fasta2), file(itr_fasta), file(contigs) from palmem_ch

	output:
	tuple val(pair_id), file("${pair_id}_discord_1.fasta"), file("${pair_id}_discord_2.fasta"), file(contigs) into discord_reads_ch

	"""
  get_discordant_reads.py ${fasta1} ${fasta2} ${itr_fasta} ${pair_id}

	rm ${fasta1}
	rm ${fasta2}
	"""
}

/*
 * Step 5: Map discordant reads to assembly
 */
process mapReads {
  label 'map_reads'

	input:
	tuple val(pair_id), file(discord1), file(discord2), file(contigs) from discord_reads_ch

	output:
	tuple val(pair_id), file(contigs), file("${pair_id}.sam") into samfile_ch

	"""
  mkdir -p bowtie2_db
	bowtie2-build ${contigs} bowtie2_db/${pair_id}_contigs
	bowtie2 --very-sensitive-local -x bowtie2_db/${pair_id}_contigs -U ${discord1},${discord2} -S ${pair_id}.sam -p ${task.cpus} -f
  rm ${discord1}
  rm ${discord2}
  rm bowtie2_db/${pair_id}_contigs*
  """
}

process filterMappedReads {
  label 'filter_mapped_reads'

  input:
  tuple val(pair_id), file(contigs), file(sam_file) from samfile_ch

  output:
  tuple val(pair_id), file(contigs), file("${pair_id}.sam.mapped.sorted") into mapped_ch

  """
	samtools view -S -b ${sam_file} -@ ${task.cpus} > ${pair_id}.bam
	rm ${sam_file}
	samtools view -b -F 4 ${pair_id}.bam -@ ${task.cpus} > ${pair_id}.bam.mapped
	rm ${pair_id}.bam
	samtools sort ${pair_id}.bam.mapped -o ${pair_id}.bam.mapped.sorted -@ ${task.cpus}
	rm ${pair_id}.bam.mapped
	samtools index ${pair_id}.bam.mapped.sorted
	samtools view ${pair_id}.bam.mapped.sorted > ${pair_id}.sam.mapped.sorted
	rm ${pair_id}.bam.mapped.sorted*
	"""
}

/*
 * Step 6: Get contigs with insertion sites
 */
process getContigsWithIRs {
  label 'contig_ir'

	input:
	tuple val(pair_id), file(contigs), file(sam_file) from mapped_ch

	output:
	tuple val(pair_id), file(sam_file), file("${pair_id}_contig_sites.fasta") into contig_ir_ch

	"""
	get_contigs_with_sites.py ${contigs} ${sam_file} ${pair_id}_contig_sites.fasta
  """
}

process getProteins {
  label 'proteins'

  input:
  tuple val(pair_id), file(sam_file), file(contig_ir) from contig_ir_ch
  path db from hmm_dir

  output:
  tuple val(pair_id), file(sam_file), file(contig_ir), file("${pair_id}_transposase.out") into proteins_ch

  """
  prodigal -i ${contig_ir} -f gff -o ${pair_id}_contigs_prodigal.gff -a ${pair_id}_translations.faa -d ${pair_id}_genes.fna -p meta
  hmmsearch --tblout ${pair_id}_transposase.out -E 1e-5 --cpu ${task.cpus} $db/$hmm_name ${pair_id}_translations.faa > /dev/null
  """
}

process getISs {
  label 'get_is'

  input:
  tuple val(pair_id), file(sam_file), file(contig_ir), file(hmm_out) from proteins_ch

  output:
  file("${pair_id}_contig_is_sites.fasta") into contig_results
  file("${pair_id}_contig_is_sites.tab") into tab_results

  """
  get_transposase.py ${hmm_out} ${sam_file} ${contig_ir} ${pair_id}
  """
}

/*
 * Get result
 */
contig_results.subscribe { it ->
    log.info "Copying results: ${contig_result_path}/${it.name}"
    it.copyTo(contig_result_path)
    }

tab_results.subscribe { it ->
    log.info "Copying results: ${contig_result_path}/${it.name}"
    it.copyTo(contig_result_path)
    }
