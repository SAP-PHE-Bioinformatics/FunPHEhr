process BARRNAP_ITS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/pybarrnap:0.5.0--pyhdfd78af_0"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.gff"), emit: gff
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    //db         = dbname ? "${dbname}" : 'euk'
    input    = fasta =~ /\.gz$/ ? fasta.name.take(fasta.name.lastIndexOf('.')) : fasta
    gunzip   = fasta =~ /\.gz$/ ? "gunzip -c ${fasta} > ${input}" : ""

    """
    $gunzip

    barrnap_its.py -i ${input} -o {prefix} -cpu $task.cpus -which 'all' -name ${meta.id}

    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        barrnap: \$(echo \$(barrnap --version 2>&1) | sed 's/barrnap//; s/Using.*\$//' )
    END_VERSIONS
    """
}