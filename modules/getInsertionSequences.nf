process getInsertionSequences {

    input:
    tuple val(sample_id), path(contig_fasta), path(contig_info), path(sam_file), path(clstr_file)

    output:
    path("${output_txt}"), emit: txt
    path("${output_fasta}"), emit: fasta

    script:
    min_is_len = params.min_is_len
    max_is_len = params.max_is_len
    output_txt="${sample_id}_insertion_sequences_info.txt"
    output_fasta="${sample_id}_insertion_sequences.fasta"
    
    """
    get_insertion_sequences.py \
        --contig_fasta ${contig_fasta} \
        --contig_info ${contig_info} \
        --sam_file ${sam_file} \
        --clstr_file ${clstr_file} \
        --min_is_len ${min_is_len} \
        --max_is_len ${max_is_len} \
        --output_prefix ${sample_id}
    """
}
