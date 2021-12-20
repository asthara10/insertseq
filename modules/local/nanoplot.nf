// Import generic module functions
include { saveFiles; getProcessName; getSoftwareName } from './functions'

params.options = [:]

process NANOPLOT {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'nanoplot', meta:meta, publish_by_meta:['id']) }
    
    conda (params.enable_conda ? "conda-forge::python=1.39.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/nanoplot:1.39.0--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/nanoplot:1.39.0--pyhdfd78af_0"
    }

    input:
    tuple val(meta), path(ontfile)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.png") , emit: png
    tuple val(meta), path("*.txt") , emit: txt
    tuple val(meta), path("*.log") , emit: log
    path  "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def input_file = ("$ontfile".endsWith(".fastq.gz") | 
        "$ontfile".endsWith(".fastq") | 
        "$ontfile".endsWith(".fq.gz") | 
        "$ontfile".endsWith(".fq")) ? "--fastq ${ontfile}" :
        ("$ontfile".endsWith(".txt")) ? "--summary ${ontfile}" : ''
    """
    NanoPlot \\
        $args \\
        -t $task.cpus \\
        $input_file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: \$(echo \$(NanoPlot --version 2>&1) | sed 's/^.*NanoPlot //; s/ .*\$//')
    END_VERSIONS
    """
}
