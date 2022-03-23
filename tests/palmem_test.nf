nextflow.enable.dsl=2

// import modules
include { palmem } from '../modules/palmem.nf'

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

workflow {
    sample_id
    .combine(fasta1)
    .combine(fasta2)
    .set { fasta_ch }

    palmem(fasta_ch)

    palmem.out.ir_1_ch
    .subscribe { it ->
        it[1].moveTo(file(params.output_dir))
        it[2].moveTo(file(params.output_dir))
    }
    palmem.out.ir_2_ch
    .subscribe { it ->
        it[1].moveTo(file(params.output_dir))
        it[2].moveTo(file(params.output_dir))
    }
    palmem.out.tab_ch
    .subscribe { it ->
        it[1].moveTo(file(params.output_dir))
    }
}
