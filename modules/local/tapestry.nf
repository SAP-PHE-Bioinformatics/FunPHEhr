process TAPESTRY {
    tag "$meta.id"
    label "process_high"
    conda 'bioconda::itsx 1.1.3'
    container "${'docker://quay.io/biocontainers/tapestry:1.0.1--pyhdfd78af_0'}"


    // The process takes a fasta file as input
    input:
    tuple val(meta), file(fasta) file(fastq)
    $telomere_seq
    // The process outputs a fasta file
    output:
    tuple val(meta), path('*_report.html')         , emit: tapestry_report, optional: true
    path "versions.yml"                             , emit: versions
    tuple val(meta), path("*.*")               , emit: results
    
    // The command to run
    script:
    def args = ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    // if fasta is gzipped, unzip it and save it to a temp file
    """
    weave -a ${fasta} -r ${fastq} -o ${prefix} -c $task.cpus -l 1000 -t ${telomere_seq}

     cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Tapestry: \$( weave -v 2>&1 | sed 's/Tapestry //g' )
    END_VERSIONS
    """
}