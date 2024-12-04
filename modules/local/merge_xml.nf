process MERGE_XML {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/pybarrnap:0.5.0--pyhdfd78af_0"

    input:
    tuple val(meta), path(xmls)

    output:
    tuple val(meta), path("${meta.id}.xml"), emit: xml

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    xmlmerge.py ${xmls} ${prefix}.xml
    """

}
