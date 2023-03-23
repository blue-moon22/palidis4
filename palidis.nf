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
include { contigCandidates } from './modules/contigCandidates.nf'
include { palmem } from './modules/palmem.nf'
include { clipIRs } from './modules/clipIRs.nf'
include { mapIRs } from './modules/mapirs.nf'
include { getInsertionSequences } from './modules/getInsertionSequences.nf'

workflow palidis {
    take:
    read_pair_ch
    contig_file_ch

    main:

    /*
     * Get maximal exact matches from reads
     */
    convertToFasta(read_pair_ch)

    palmem(convertToFasta.out)

    clipIRs(palmem.out.ir_ch)

    /*
     * Filter contigs
     */
    filterContigs(contig_file_ch)

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

    runProdigal(filterContigs.out)

    runProdigal.out
    .combine(interproscan_ch)
    .set { proteins_ch }

    runInterproscan(proteins_ch)

    contigCandidates(runInterproscan.out.fasta)

    buildDB(contigCandidates.out.ref)

    clipIRs.out
    .join(buildDB.out.contig_db_ch)
    .set { irs_contig_ch }

    /*
     * Map IRs to contigs
     */
     mapIRs(irs_contig_ch)

    /*
     * Get insertion sequences and corresponding information
     */
     contigCandidates.out.fasta_info
     .join(mapIRs.out)
     .set { contig_info_ch }

    getInsertionSequences(contig_info_ch)
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
    .map { row -> tuple(row.sample_id, file(row.read1), file(row.read2)) }
    .groupTuple()
    .set { read_pair_ch }

    Channel
    .fromPath(params.manifest, checkIfExists: true)
    .splitCsv(header:true, sep:"\t")
    .map { row -> tuple(row.sample_id, file(row.contigs_path)) }
    .set { contig_file_ch }

    palidis(read_pair_ch, contig_file_ch)

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
