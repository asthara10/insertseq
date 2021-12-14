#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/insertseq
========================================================================================
    Github : https://github.com/nf-core/insertseq
    Website: https://nf-co.re/insertseq
    Slack  : https://nfcore.slack.com/channels/insertseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { INSERTSEQ } from './workflows/insertseq'

//
// WORKFLOW: Run main nf-core/insertseq analysis pipeline
//
workflow NFCORE_INSERTSEQ {
    INSERTSEQ ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_INSERTSEQ ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
