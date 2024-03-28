

process COMPLEASM {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::0.2.4"
    container 'biocontainers/compleasm:0.2.4--pyh7cba7a3_0'

    input:
    tuple val(meta), path('tmp_input/*')
    val mode                              // Required:    One of genome, proteins, or transcriptome
    val lineage                           // Required:    lineage to check against, "auto" enables --auto-lineage instead
    path busco_lineages_path              // Recommended: path to busco lineages - downloads if not set

    output:
    tuple val(meta), path("*summary.txt")                             , emit: batch_summary
    tuple val(meta), path("short_summary.*.txt")                      , emit: short_summaries_txt, optional: true
    tuple val(meta), path("short_summary.*.json")                     , emit: short_summaries_json, optional: true
    tuple val(meta), path("*-busco/*/run_*/full_table.tsv")           , emit: full_table, optional: true
    tuple val(meta), path("*-busco/*/run_*/missing_busco_list.tsv")   , emit: missing_busco_list, optional: true
    tuple val(meta), path("*-busco/*/run_*/single_copy_proteins.faa") , emit: single_copy_proteins, optional: true
    tuple val(meta), path("*-busco/*/run_*/busco_sequences")          , emit: seq_dir
    tuple val(meta), path("*-busco/*/translated_proteins")            , emit: translated_dir, optional: true
    tuple val(meta), path("*-busco")                                  , emit: busco_dir
    path "versions.yml"                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}-${lineage}"
    def busco_lineage = lineage.equals('auto') ? '--auto-lineage' : "-l ${lineage}"
    def busco_lineage_dir = busco_lineages_path ? "--download_path ${busco_lineages_path}" : ''
    """"
    # Ensure the input is uncompressed
    INPUT_SEQS=input_seqs
    mkdir "\$INPUT_SEQS"
    cd "\$INPUT_SEQS"
    for FASTA in ../tmp_input/*; do
        if [ "\${FASTA##*.}" == 'gz' ]; then
            gzip -cdf "\$FASTA" > \$( basename "\$FASTA" .gz )
        else
            ln -s "\$FASTA" .
        fi
    done
    cd ..

    python compleasm.py run \\
        -t $task.cpus \\
        -a "\$INPUT_SEQS" \\
        -o ${prefix} \\
        --mode $mode \\
        $busco_lineage \\
        -L $busco_lineage_dir \\
        $args

    # clean up
    rm -rf "\$INPUT_SEQS"

    # Move files to avoid staging/publishing issues
    mv ${prefix}-busco/batch_summary.txt ${prefix}-busco.batch_summary.txt
    mv ${prefix}-busco/*/short_summary.*.{json,txt} . || echo "Short summaries were not available: No genes were found."

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        compleasm: \$( 0.2.4 )
    END_VERSIONS
    """
}