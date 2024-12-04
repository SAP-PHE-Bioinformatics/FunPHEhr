process FUNANNOTATE_IPRSCAN {
    tag "$meta.id"
    label 'process_medium'
    conda 'bioconda:: funannotate 1.8.15'
    container "${'biocontainers/funannotate:1.8.15--pyhdfd78af_2'}"


    // The process takes a fasta file as input
    input:
    tuple val(meta), path(files)
    // The process outputs a fasta file
    output:
    path "versions.yml"                        , emit: versions
    tuple val(meta), path("*.*")               , emit: results
    tuple val(meta), path("results/*.tsv")             , emit: tsv
    
    // The command to run
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    funannotate iprscan -i $files -m local --cpus 12 -o results

    
     cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         funannotate: \$( funannotate --version 2>&1 | sed 's/funannotate //g' )
    END_VERSIONS
    """
}