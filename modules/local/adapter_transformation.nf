// Import generic module functions
include { saveFiles; getProcessName; getSoftwareName } from './functions'

params.options = [:]

process ADAPTER_TRANSFORMATION {
    label "process_low"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    input:
    tuple val(meta), val(adapter_1F), val(payload)
    val adapter_2

    output:
    tuple val(meta), env(adapter_1m)                 , emit: adapter1
    tuple val(meta), env(adapter_2m), env(payload_rc), emit: adapter2
    path "versions.yml"                              , emit: versions

    script: 
    """
    payload_rc=`echo $payload | tr ACGTacgt TGCATGCA | sed -r 'G;:a;s/^(.)(.*\\n)/\\2\\1/;ta;s/\\n//'`
    adapter_1m=`echo $adapter_1F | tr acgt ACGT`
    adapter_2m=`echo $adapter_2 | tr acgt ACGT`

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
