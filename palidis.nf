/*
 * Nextflow pipeline for identifying insertion sequences from metagenomic data
 *
 * Author:
 * Victoria Carr vc11@sanger.ac.uk
 *
 */

nextflow.enable.dsl=2

// import modules
include { convertToFasta } from './modules/convertToFasta.nf'
include { filterContigs } from './modules/filterContigs.nf'
include { buildDB } from './modules/buildDB.nf'
include { runProdigal } from './modules/runProdigal.nf'
include { installInterproscan } from './modules/installInterproscan.nf'
include { runInterproscan } from './modules/runInterproscan.nf'
include { mergeTSV } from './modules/merge.nf'
include { palmem } from './modules/palmem.nf'
include { clipIRs } from './modules/clipIRs.nf'
include { mapIRs } from './modules/mapirs.nf'
include { getInsertionSequences } from './modules/getInsertionSequences.nf'
include { getCandidateInsertionSequences } from './modules/getCandidateInsertionSequences.nf'

workflow palidis {
    take:
    reads_ch
    contig_file_ch

    main:

    /*
     * Get maximal exact matches from reads
     */
    convertToFasta(reads_ch)

    palmem(convertToFasta.out)

    clipIRs(palmem.out)

    /*
     * Filter contigs
     */
    filterContigs(contig_file_ch)

    buildDB(filterContigs.out.db_ch)

    /*
     * Map IRs to candidate contigs
     */
    clipIRs.out
    .join(buildDB.out)
    .set { irs_contig_ch }

    mapIRs(irs_contig_ch)

    /*
     * Get candidate Insertion Sequences
     */
    filterContigs.out.fasta_ch
    .join(mapIRs.out)
    .set { candidate_contigs_ch }

    getCandidateInsertionSequences(candidate_contigs_ch)

    /*
     * Download Interproscan
     */
    if (!file("${params.db_path}/${params.interproscan_db}").exists()) {
        installInterproscan()

        db_path = file("${params.db_path}")
        db_path.mkdir()

        installInterproscan.out
        .set { interproscan_ch }

        interproscan_ch
        .subscribe{ it ->
            it.copyTo("${db_path}")
        }

    } else {
        Channel
        .fromPath(file("${params.db_path}/${params.interproscan_db}"))
        .set { interproscan_ch }
    }

    /*
     * Annotate transposase
     */
    runProdigal(getCandidateInsertionSequences.out.fasta_ch)

    runProdigal.out.prot
    .splitFasta(by: params.chunk_size, file:true)
    .combine(interproscan_ch)
    .set { chunk_ch }

    runInterproscan(chunk_ch)

    mergeTSV(runInterproscan.out.groupTuple(size: 0))

    getCandidateInsertionSequences.out.fasta_txt_ch
    .join(mergeTSV.out)
    .set { contig_tsv_ch }

    /*
     * Get insertion sequences and corresponding information
     */
    getInsertionSequences(contig_tsv_ch)
    is_info_ch = getInsertionSequences.out.txt
    is_fasta_ch = getInsertionSequences.out.fasta

    emit:
    is_fasta_ch
    is_info_ch
}

workflow {
    // Define parameters
    batch_path = file("./${params.batch_name}")
    batch_path.mkdir()

    /*
     * Parameters
     */
    Channel
    .fromPath(params.manifest, checkIfExists: true)
    .splitCsv(header:true, sep:"\t")
    .map { row -> tuple(row.sample_id, file(row.read_directory)) }
    .groupTuple()
    .set { reads_ch }

    Channel
    .fromPath(params.manifest, checkIfExists: true)
    .splitCsv(header:true, sep:"\t")
    .map { row -> tuple(row.sample_id, file(row.contigs_path)) }
    .set { contig_file_ch }

    palidis(reads_ch, contig_file_ch)

    // Publish IS fasta sequences
    palidis.out.is_fasta_ch
    .subscribe { it ->
        it.copyTo("${batch_path}")
    }

    // Publish annotations
    palidis.out.is_info_ch
    .subscribe { it ->
        it.copyTo("${batch_path}")
    }
}
