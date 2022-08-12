/*
 * Build bowtie2 database
 */
process buildDB {

    input:
    tuple val(sample_id), path(contigs_path)

    output:
    tuple val(sample_id), path("contigs_db.tar"), emit: contig_db1_ch
    tuple val(sample_id), path("contigs_db.tar"), emit: contig_db2_ch

    """
    bowtie2-build ${contigs_path} ${sample_id}_contigs
    tar -cf contigs_db.tar ${sample_id}_contigs*
    """
}
