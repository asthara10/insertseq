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
    sequences                                 // channel: [ val(meta), [ adapter_1F, payload ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ adapter_1F, payload ] ]
def create_fastq_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.fastq

    def array = []
    if (!file(row.adapter_1F).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Adapter_1F sequence does not exist!\n${row.adapter_1F}"
    }
    if (!file(row.payload).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Payload primer sequence does not exist!\n${row.payload}"
    }
    array = [ meta, [ file(row.adapter_1F), file(row.payload) ] ]

    return array
}
