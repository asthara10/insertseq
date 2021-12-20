//
// Check input samplesheet and get read channels
//

params.options = [:]

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check' addParams( options: params.options )

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channels(it) }
        .set { sequences }

    emit:
    sequences                                 // channel: [ val(meta), [ fastq, adapter_1F, payload ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ adapter_1F, payload ] ]
def create_fastq_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id  = row.fastq.split('/')[-1].split('.f')[0]

    def array = []
    if (!file(row.fastq).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Input file does not exist!\n${row.fastq}"
    }
    array = [ meta, [ file(row.fastq), row.adapter_1F, row.payload ] ]

    return array
}
