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
include { palmem } from './modules/palmem.nf'
include { mapReads as mapReads1 } from './modules/mapreads.nf'
include { mapReads as mapReads2 } from './modules/mapreads.nf'
include { getCandidateITRs } from './modules/getCandidateITRs.nf'
include { clusterReads } from './modules/clusterReads.nf'
include { getITRs } from './modules/getITRs.nf'
include { buildBLASTDB } from './modules/buildBLASTDB.nf'
include { searchISfinder } from './modules/searchISfinder.nf'
include { searchCOBSIndex } from './modules/searchCOBSIndex.nf'
include { getISInfoWithCOBS; getISInfoWithoutCOBS } from './modules/getISInfo.nf'

workflow get_IS_annotations {
    take:
    read_pair_ch
    contig_file_ch

    main:
    convertToFasta(read_pair_ch)

    palmem(convertToFasta.out)

    filterContigs(contig_file_ch)

    buildDB(filterContigs.out)

    /*
     * Map reads to contigs
     */
    Channel
    .of('1')
    .set { pair1_ch }

    palmem.out.ir_1_ch
    .join(buildDB.out.contig_db1_ch)
    .combine(pair1_ch)
    .set { contigs_reads1_ch }

    mapReads1(contigs_reads1_ch)
    mapping1_ch = mapReads1.out

    Channel
    .of('2')
    .set { pair2_ch }

    palmem.out.ir_2_ch
    .join(buildDB.out.contig_db2_ch)
    .combine(pair2_ch)
    .set { contigs_reads2_ch }

    mapReads2(contigs_reads2_ch)
    mapping2_ch = mapReads2.out

    /*
     * Get contigs and reads with candidate ITRs
     */
     contig_file_ch
     .join(mapping1_ch)
     .join(mapping2_ch)
     .join(convertToFasta.out)
     .join(palmem.out.tab_ch)
     .set { mapping_contigs_ch }

    getCandidateITRs(mapping_contigs_ch)

    reads_itrs_ch = getCandidateITRs.out.reads_itrs_ch
    tab_ch = getCandidateITRs.out.tab_ch

    clusterReads(reads_itrs_ch)
    cluster_ch = clusterReads.out.cluster_ch

    cluster_ch
    .join(tab_ch)
    .join(filterContigs.out)
    .set { into_get_itr_ch }

    getITRs(into_get_itr_ch)
    is_tab_ch = getITRs.out.is_tab_ch
    is_fasta_ch1 = getITRs.out.is_fasta_ch1
    is_fasta_ch2 = getITRs.out.is_fasta_ch2
    is_fasta_ch = getITRs.out.is_fasta_ch

    Channel
    .fromPath(params.isfinder_seq, checkIfExists: true)
    .set { fasta_db_ch }
    buildBLASTDB(fasta_db_ch)
    blast_db_ch = buildBLASTDB.out

    is_fasta_ch1
    .combine(blast_db_ch)
    .set{ is_seq_ch }
    searchISfinder(is_seq_ch)
    blast_out_ch = searchISfinder.out

    Channel
    .fromPath(params.isfinder_info_csv, checkIfExists: true)
    .set { isfinder_info_ch }

    if (params.cobs_index) {
        Channel
        .fromPath(params.cobs_index, checkIfExists: true)
        .set { cobs_index_ch }

        is_fasta_ch2
        .combine(cobs_index_ch)
        .set{ cobs_seq_ch }

        searchCOBSIndex(cobs_seq_ch)
        cobs_out_ch = searchCOBSIndex.out

        is_tab_ch
        .join(blast_out_ch)
        .join(cobs_out_ch)
        .set { is_annot_ch }

        getISInfoWithCOBS(is_annot_ch, isfinder_info_ch)
        is_info_ch = getISInfoWithCOBS.out
    } else {
        is_tab_ch
        .join(blast_out_ch)
        .set { is_annot_ch }

        getISInfoWithoutCOBS(is_annot_ch, isfinder_info_ch)
        is_info_ch = getISInfoWithoutCOBS.out
    }

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

    get_IS_annotations(read_pair_ch, contig_file_ch)

    // Publish IS fasta sequences
    get_IS_annotations.out.is_fasta_ch
    .subscribe { it ->
        it.copyTo("${batch_path}")
    }

    // Publish annotations
    get_IS_annotations.out.is_info_ch
    .subscribe { it ->
        it.copyTo("${batch_path}")
    }
}
