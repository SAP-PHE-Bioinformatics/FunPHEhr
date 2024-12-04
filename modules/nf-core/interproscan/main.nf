process INTERPROSCAN {
    tag "$meta.id"
    label 'process_medium_mem'

    conda "${moduleDir}/environment.yml"
    container "docker://interpro/interproscan:latest"
    clusterOptions = " --nodelist=frgeneseq-control"

    input:
    tuple val(meta), path(fasta)
    path(interproscan_database, stageAs: 'data')

    output:
    tuple val(meta), path('*.tsv') , optional: true, emit: tsv
    tuple val(meta), path('*.xml') , optional: true, emit: xml
    tuple val(meta), path('*.gff3'), optional: true, emit: gff3
    tuple val(meta), path('*.json'), optional: true, emit: json
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.name.endsWith(".gz")
    def fasta_name = fasta.name.replace(".gz", "")
    def file_name = fasta_name.replace(".fasta", "")
    """
    if [ -d 'data' ]; then
        # Find interproscan.properties to link data/ from work directory
        INTERPROSCAN_DIR="./data"
        INTERPROSCAN_PROPERTIES="./data/interproscan.properties"
        cp "\$INTERPROSCAN_PROPERTIES" .
        sed -i "/^bin\\.directory=/ s|.*|bin.directory=\$INTERPROSCAN_DIR/bin|" interproscan.properties
        export INTERPROSCAN_CONF=interproscan.properties
    fi # else use sample DB included with conda ( testing only! )

    if ${is_compressed} ; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    ./data/interproscan.sh \\
        --cpu ${task.cpus} \\
        --input ${fasta_name} \\
        ${args} \\
        --output-file-base ${file_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        interproscan: \$( interproscan.sh --version | sed '1!d; s/.*version //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.{tsv,xml,json,gff3}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        interproscan: \$( interproscan.sh --version | sed '1!d; s/.*version //' )
    END_VERSIONS
    """
}