// Import generic module functions
include { saveFiles; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ADAPTER_TRANSFORMATION {
    label "process_low"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }

    input:
    tuple val(meta), val(adapter_1F), val(payload)
    params.adapter_2

    output:
    tuple val(meta), env(adapter_1m),                , emit: adapter1
    tuple val(meta), env(adapter_2m), env(payload_rc), emit: adapter2
    path "versions.yml",                               emit: versions

    script: 
    """
    payload_rc=`echo $payload | tr ACGTacgt TGCATGCA | rev`
    adapter_1m=`echo $adapter_1F | tr acgt ACGT`
    adapter_2m=`echo $adapter_2 | tr acgt ACGT`

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
    END_VERSIONS
    """
}
