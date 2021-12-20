/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def valid_params = [
    umi_type       : ['TV', 'NYR'],
 ]

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowInsertseq.initialise(params, log, valid_params)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.genome, params.ref_vector, params.repeats ]

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Loaded from modules/local/
//
// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

include{ ADAPTER_TRANSFORMATION } from '../modules/local/adapter_transformation' addParams( options: [:] )
include { NANOPLOT } from '../modules/local/nanoplot'  addParams( options: [:] )
include { NANOFILT } from '../modules/local/nanofilt'  addParams( options: [:] )

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check' addParams( options: [:] )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//def multiqc_options   = modules['multiqc']
//multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// MODULE: Installed directly from nf-core/modules
//



include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'  addParams( options: [publish_files : ['_versions.yml':'']] )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow INSERTSEQ {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    .sequences
    .multiMap {
        meta, value ->
            all_values = value.flatten()
            files: [meta, all_values[0]]
            adapters: [meta, all_values[1], all_values[2]]
    }
    .set { ch_inputs }
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    ch_inputs.files.multiMap { it -> plot: filter: it }.set { input }

    //
    // MODULE: prepare adapter sequences
    //
    ADAPTER_TRANSFORMATION (
        ch_inputs.adapters,
        params.adapter_2
    )
    adapter1 = ADAPTER_TRANSFORMATION.out.adapter1
    adapter2 = ADAPTER_TRANSFORMATION.out.adapter2
    ch_versions = ch_versions.mix(ADAPTER_TRANSFORMATION.out.versions.first().ifEmpty(null))

    //
    // MODULE: run NanoPlot to assess read quality and length
    //
    NANOPLOT {
        input.plot
    }
    .set { ch_nanoplot }
    ch_versions = ch_versions.mix(NANOPLOT.out.versions.first().ifEmpty(null))

    //
    // MODULE: filter reads by quality
    //
    NANOFILT {
        input.filter
    }


    /*CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowInsertseq.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    */
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
