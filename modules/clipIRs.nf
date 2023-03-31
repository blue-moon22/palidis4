/*
 * Clip reads and cluster inverted repeats
 */
process clipIRs {

    input:
    tuple val(sample_id), file(ir_fasta), file(tab_file)

    output:
    tuple val(sample_id), path(output_fasta)

    script:
    output_fasta="${sample_id}_clipped_irs.fasta"

    """
    clip_reads.py --read_fasta ${ir_fasta} --tab_file ${tab_file} --output_prefix ${sample_id}_clipped
    """
}
