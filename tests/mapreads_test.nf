nextflow.enable.dsl=2

// import modules
include { mapReads } from '../modules/mapreads.nf'

// Input files
Channel
.value(params.sample_id)
.set{ sample_id }

Channel
.fromPath(file(params.fasta1), checkIfExists:true)
.set{ fasta1 }

Channel
.fromPath(file(params.fasta2), checkIfExists:true)
.set{ fasta2 }

Channel
.fromPath(file(params.db_path), checkIfExists:true)
.set { db_path }

Channel
.of('1')
.set { pair }

workflow {
    sample_id
    .combine(fasta1)
    .combine(fasta2)
    .combine(db_path)
    .combine(pair)
    .set { contigs_reads1_ch }

    mapReads(contigs_reads1_ch)

    mapReads.out
    .subscribe { it ->
        it[2].moveTo(file(params.output_dir))
    }
}
