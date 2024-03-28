process FUNANNOTATE_MASK {
    tag "$meta.id"
    label 'process_medium'
    conda 'bioconda:: funannotate 1.8.15'
    container "${'docker://quay.io/biocontainers/funannotate:1.8.15--pyhdfd78af_2'}"


    // The process takes a fasta file as input
    input:
    tuple val(meta), file(fasta)
    // The process outputs a fasta file
    output:
    tuple val(meta), path('*_masked.fasta')         , emit: masked
    path "versions.yml"                             , emit: versions
    
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
    funannotate mask -i ${prefix}.fasta -o ${prefix}_masked.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         funannotate: \$( funannotate --version 2>&1 | sed 's/funannotate //g' )
    END_VERSIONS
    """
}