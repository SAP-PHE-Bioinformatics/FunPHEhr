// itsx nf-core module

// The docker image used

process ITSX {
    label "ITSx"
    conda 'bioconda::itsx 1.1.3'
    container "${'docker://quay.io/biocontainers/itsx:1.1.3--py27_0'}"


    // The process takes a fasta file as input
    input:
    tuple val(meta), file(fasta)
    // The process outputs a fasta file
    output:
    tuple val(meta), path('*_.full.fasta')         , emit: itsx_full, optional: true
    tuple val(meta), path('*_.summary.txt')         , emit: istx_summary, optional: true
    path "versions.yml"                             , emit: versions
    
    // The command to run
    script:
    def args = ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    ITSx -i ${fasta} -o ${prefix} --save_regions all --cpu 8
    """
}