process EXTRACTKRAKENTOPMATCH {

    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krakentools:1.2--pyh5e36f6f_0':
        'biocontainers/krakentools:1.2--pyh5e36f6f_0'}"

    input:
    tuple val(meta), path(report), path(classified_reads_fastq), path(classified_reads_assignment) 

    output:
    tuple val(meta), path("*.{fastq,fasta}"), emit: extracted_kraken2_reads

    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--fastq-output") ? "fastq" : "fasta"
    def input_reads_command = meta.single_end ? "-s $classified_reads_fastq" : "-s1 ${classified_reads_fastq[0]} -s2 ${classified_reads_fastq[1]}"
    def output_reads_command = meta.single_end ? "-o ${prefix}.extracted_kraken2_read.${extension}" : "-o ${prefix}.extracted_kraken2_read_1.${extension} -o2 ${prefix}.extracted_kraken2_read_2.${extension}"
    def report_option = report ? "-r ${report}" : ""
    def VERSION = '1.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    
    extract_kraken_top_match.py \\
        ${args} \\
        -k $classified_reads_assignment \\
        $report_option \\
        $input_reads_command \\
        $output_reads_command \\
        --append

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        extract_kraken_reads.py: ${VERSION}
    END_VERSIONS
    """
}