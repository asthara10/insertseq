// Import generic module functions
include { saveFiles; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ADAPTER_TRANSFORMATION {
    tag "$samplesheet"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    input:
    path samplesheet

    output:
    path '*.csv'       , emit: csv
    path "versions.yml", emit: versions

    script: // This script is bundled with the pipeline, in nf-core/insertseq/bin/
    """
    ontarget_rc=`echo $ontarget | tr ACGTacgt TGCATGCA | rev`
    adapter_1m=`echo $adapter_1F | tr acgt ACGT`
    adapter_2m=`echo $adapter_2 | tr acgt ACGT`
    
    check_samplesheet.py \\
        $samplesheet \\
        samplesheet.valid.csv

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
