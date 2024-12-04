// itsx nf-core module

// The docker image used

process AURICLASS {
    tag "$meta.id"
    label 'process_medium'
    conda 'bioconda::itsx 1.1.3'
    container "${'docker://quay.io/biocontainers/auriclass:0.5.4--pyhdfd78af_0'}"


    // The process takes a fasta file as input
    input:
    tuple val(meta), file(fasta)
    // The process outputs a fasta file
    output:
    tuple val(meta), path('*.tsv')         , emit: auriclass_report
    path "versions.yml"                             , emit: versions
    
    // The command to run
    script:
    def args = ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    auriclass -o ${prefix}.cauris_clade.tsv ${fasta}


     cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        AuriClass : \$( auriclass --version )
    END_VERSIONS
    """
}