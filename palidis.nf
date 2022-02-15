/*
 * Nextflow pipeline for identifying insertion sequences from metagenomic data
 *
 * Author:
 * Victoria Carr <victoriacarr018@gmail.com
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
include { createITRCatalog } from './modules/createITRCatalog.nf'
include { assignITRClusters } from './modules/assignITRClusters.nf'
include { createISCatalog } from './modules/createISCatalog.nf'

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
    itr_fasta_ch = clusterReads.out.clipped_read_ch

    cluster_ch
    .join(tab_ch)
    .join(filterContigs.out)
    .set { into_get_itr_ch }

    getITRs(into_get_itr_ch)
    itr_clusters_ch = getITRs.out.itr_clusters_ch
    is_tab_ch = getITRs.out.is_tab_ch

    emit:
    itr_fasta_ch
    itr_clusters_ch
    is_tab_ch
}

workflow create_IS_catalog {
    take:
    itrs_fasta_ch
    is_annot_ch
    itr_clusters_ch

    main:
    createITRCatalog(itrs_fasta_ch)
    all_itrs_ch = createITRCatalog.out

    assignITRClusters(all_itrs_ch)
    all_itrs_tab_ch = assignITRClusters.out

    createISCatalog(all_itrs_tab_ch, itr_clusters_ch, is_annot_ch)
    is_catalog_ch = createISCatalog.out

    emit:
    is_catalog_ch
}

workflow {
    // Define parameters
    batch_path = file("./${params.batch_name}")
    batch_path.mkdir()

    if (params.get_IS_annotations) {
        /*
         * Parameters
         */
        Channel
        .fromPath(params.manifest)
        .splitCsv(header:true, sep:"\t")
        .map { row -> tuple(row.sample_id, file(row.read1), file(row.read2)) }
        .groupTuple()
        .set { read_pair_ch }

        Channel
        .fromPath(params.manifest)
        .splitCsv(header:true, sep:"\t")
        .map { row -> tuple(row.sample_id, file(row.contigs_path)) }
        .set { contig_file_ch }

        get_IS_annotations(read_pair_ch, contig_file_ch)

        // Publish batch of candidate ITRs
        get_IS_annotations.out.itr_fasta_ch
        .flatten()
        .subscribe { it ->
            it.copyTo("${batch_path}")
        }

        // Publish itr clusters file for batch
        get_IS_annotations.out.itr_clusters_ch
        .subscribe { it ->
            it.copyTo("${batch_path}")
        }

        // Publish tab file for batch
        get_IS_annotations.out.is_tab_ch
        .subscribe { it ->
            it.copyTo("${batch_path}")
        }
    }

    if (params.create_catalog) {

        Channel
        .fromPath("${batch_path}/*_ITRs.fasta", checkIfExists:true)
        .collect()
        .set { itrs_fasta_ch }

        Channel
        .fromPath("${batch_path}/*_insertion_sequence_annotations.tab", checkIfExists:true)
        .collect()
        .set { is_annot_ch }

        Channel
        .fromPath("${batch_path}/*_reads_itr_clusters.txt", checkIfExists:true)
        .collect()
        .set { itr_clusters_ch }

        create_IS_catalog(itrs_fasta_ch, is_annot_ch, itr_clusters_ch)

        // Publish catalog
        create_IS_catalog.out.is_catalog_ch
        .subscribe { it ->
            it.copyTo("./")
        }
    }
}
