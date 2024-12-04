// itsx nf-core module

// The docker image used

process ITSX {
    tag "$meta.id"
    label 'process_medium'
    conda 'bioconda::itsx 1.1.3'
    container "${'docker://quay.io/biocontainers/itsx:1.1b--1'}"


    // The process takes a fasta file as input
    input:
    tuple val(meta), file(fasta)
    // The process outputs a fasta file
    output:
    tuple val(meta), path('*.full.fasta')         , emit: itsx_full
    tuple val(meta), path('*.summary.txt')         , emit: istx_summary, optional: true
    path "versions.yml"                             , emit: versions
    tuple val(meta), path("*.*")               , emit: results
    
    // The command to run
    script:
    def args = ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    // if fasta is gzipped, unzip it
    """
    if [[ ${fasta} == *.gz ]]; then
        gzip -cdf ${fasta} > ${prefix}.fasta
    else
        cp ${fasta} ${prefix}.fasta
    fi
    ITSx -i ${prefix}.fasta -o ${prefix} --save_regions all --cpu $task.cpus 


     cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ITSx: 1.1.3
    END_VERSIONS
    """
}