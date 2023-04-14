process getCandidateInsertionSequences {

    input:
    tuple val(sample_id), path(contig_fasta), path(sam_file)

    output:
    tuple val(sample_id), path(output_fasta), emit: fasta_ch
    tuple val(sample_id), path(output_fasta), path(output_txt), emit: fasta_txt_ch

    script:
    min_is_len = params.min_is_len
    max_is_len = params.max_is_len
    output_txt="${sample_id}_candidate_insertion_sequences_info.txt"
    output_fasta="${sample_id}_candidate_insertion_sequences.fasta"

    """
    get_candidate_insertion_sequences.py \
        --contig_fasta ${contig_fasta} \
        --sam_file ${sam_file} \
        --min_is_len ${min_is_len} \
        --max_is_len ${max_is_len} \
        --output_prefix ${sample_id}
    """
}
