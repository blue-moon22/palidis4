/*
 * Clip reads and cluster inverted repeats
 */
process clipIRs {

    input:
    tuple val(sample_id), file(ir_1), file(ir_2), file(tab_file)

    output:
    tuple val(sample_id), path(output_fasta)

    script:
    output_fasta="${sample_id}_clipped_irs.fasta"

    """
    cat ${ir_1} ${ir_2} > ${sample_id}_IR.fasta
    clip_reads.py --read_fasta ${sample_id}_IR.fasta --tab_file ${tab_file} --output_prefix ${sample_id}_clipped
    """
}
