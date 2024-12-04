process HELIXER {
    // use conatiner in singulatirty from docker://gglyptodon/helixer-docker:helixer_v0.3.2_cuda_11.8.0-cudnn8
    //tool github https://github.com/weberlab-hhu/Helixer
    tag "$meta.id"
    label 'gpu_medium'

    conda 'helixer=0.3.2'
    container 'docker://gglyptodon/helixer-docker:helixer_v0.3.2_cuda_11.8.0-cudnn8'

    input:
    tuple val(meta), file(fasta)
    path(model)

    output:
    tuple val(meta), path('*_helixer.gff3')         , emit: gff3
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                    = task.ext.args ?: ''
    def prefix                  = task.ext.prefix ?: "${meta.id}"
    def species                 = task.ext.species ?: "${meta.species}"
    //replace any space in species with _
    species                     = species.replaceAll("\\s", "_")


    """
    Helixer.py \
        --model-filepath $model \
        --subsequence-length 21384 \
        --fasta-path ${ fasta } \
        --compression gzip \
        --species ${species} \
        --gff-output-path ${prefix}_helixer.gff3 \

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Helixer: \$( Helixer.py --version 2>&1 | sed 's/Helixer //g' )
    END_VERSIONS
    """



    
}

