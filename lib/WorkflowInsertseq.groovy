//
// This file holds several functions specific to the workflow/insertseq.nf in the nf-core/insertseq pipeline
//

class WorkflowInsertseq {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log, valid_params) {
        genomeExistsError(params, log)

        // Genome
        if (!params.genome) {
            log.error "Genome fasta file not specified with e.g. '--genome genome.fa' or via a detectable config file."
            System.exit(1)
        }

        // Adapters
        if (!params.adapter_1F) {
            log.error "First adapter sequence not specified with e.g. '--adapter_1 ACTG' or via a detectable config file."
            System.exit(1)
        }
        if (!params.adapter_2) {
            log.error "Second adapter sequence not specified with e.g. '--adapter_2 ACTG' or via a detectable config file."
            System.exit(1)
        }
        if (!params.payload) {
            log.error "Payload primer sequence not specified with e.g. '--payload ACTG' or via a detectable config file."
            System.exit(1)
        }

        if (!params.skip_umiclustering) {
            if (!valid_params["umi_type"].contains(params.umi_type)) {
                log.error "Invalid option: '${params.umi_type}'. Valid options for '--umi_type': ${valid_params['umi_type'].join(', ')}."
                System.exit(1)
            }
        }
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    //
    // Exit pipeline if incorrect --genome key provided
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            log.error "=============================================================================\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "==================================================================================="
            System.exit(1)
        }
    }
}
