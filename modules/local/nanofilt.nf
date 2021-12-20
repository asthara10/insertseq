// Import generic module functions
include { saveFiles; getProcessName; getSoftwareName } from './functions'

params.options = [:]

process NANOFILT {
    label "process_low"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'filtering', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::nanofilt=2.8.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/nanofilt:2.8.0--py_0"
    } else {
        container "quay.io/biocontainers/nanofilt:2.8.0--py_0"
    }

    input:
    tuple val(meta), val(input_file)

    output:
    tuple val(meta), file("*_qualityFilt.fastq"), emit: filtered
    path "versions.yml"                         , emit: versions

    script: 
    def args = task.ext.args ?: ''
    """
    NanoFilt \\
        $args \\
        $input_file > \\
        ${meta.id}_qualityFilt.fastq

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        nanofilt: \$(NanoFilt --version | sed 's/NanoFilt //g')
    END_VERSIONS
    """
}
