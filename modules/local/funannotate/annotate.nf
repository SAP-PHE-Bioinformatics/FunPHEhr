process ANNOTATE {
    tag "$meta.id"
    label 'process_medium'
    conda 'bioconda:: funannotate 1.8.15'
    container 'docker://nextgenusfs/funannotate'


    // The process takes a fasta file as input
    input:
    tuple val(meta), path(fasta), path(gff), path(interpro), path(eggnog)
    // The process outputs a fasta file
    output:
    path "versions.yml"                        , emit: versions
    tuple val(meta), path("${meta.id}_annotate")               , emit: results
    
    // The command to run
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def species = "${meta.species}" ?: ''
    real_zip = "readlink ${eggnog}"
    """

    mkdir -p ./${prefix}_annotate
    
    funannotate annotate --fasta ${fasta} -o ./${prefix}_annotate --cpus $task.cpus -s ${species} --gff ${gff} --iprscan ${interpro} \
    --eggnog ${eggnog} --header_length 35 \
    --isolate ${prefix}



     cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         funannotate: \$( funannotate --version 2>&1 | sed 's/funannotate //g' )
    END_VERSIONS
    """
}