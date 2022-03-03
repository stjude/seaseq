version 1.0

task sicer {

    input {
        File bedfile
        File? control_bed
        File chromsizes
        String default_location = "PEAKS_files/BROAD_peaks"
        String coverage_location = "COVERAGE_files/NARROW_peaks"

        Int redundancy = 1
        Int window = 200
        Int fragment_size = 150
        Float genome_fraction = 0.86
        Int gap_size = 200
        Int evalue = 100

        String outputname = basename(bedfile,'.bed')

        Int memory_gb = 5
        Int max_retries = 1
        Int ncpu = 20
    }

    command <<<
        mkdir -p ~{default_location} ~{coverage_location}

        sicer \
            -cpu ~{ncpu} \
            -t ~{bedfile} \
            ~{if defined(control_bed) then "-c" + control_bed else ""} \
            -sc ~{chromsizes} \
            -rt ~{redundancy} \
            -w ~{window} \
            -f ~{fragment_size} \
            -egf ~{genome_fraction} \
            -g ~{gap_size} \
            -e ~{evalue}

        mv ~{basename(bedfile,'.bed')}-W~{window}-G~{gap_size}.scoreisland ~{outputname}-W~{window}-G~{gap_size}.scoreisland
        mv ~{basename(bedfile,'.bed')}-W~{window}-normalized.wig ~{outputname}-W~{window}-normalized.wig
        mv ~{basename(bedfile,'.bed')}-W~{window}-G~{gap_size}-islands-summary ~{outputname}-W~{window}-G~{gap_size}-islands-summary
        if [ -f "~{basename(bedfile,'.bed')}-W~{window}-G~{gap_size}-FDR0.01-island.bed" ]; then
            mv ~{basename(bedfile,'.bed')}-W~{window}-G~{gap_size}-FDR0.01-island.bed ~{outputname}-W~{window}-G~{gap_size}-FDR0.01-island.bed
        fi
        gzip *wig
        mv ~{outputname}-W~{window}-normalized.wig.gz ~{coverage_location}
        mv ~{outputname}-* ~{default_location}
    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'abralab/sicer:v1.1.0'
        cpu: ncpu
    }
    output {
        File scoreisland = "~{default_location}/~{outputname}-W~{window}-G~{gap_size}.scoreisland"
        File wigfile = "~{coverage_location}/~{outputname}-W~{window}-normalized.wig.gz"
        File? summary = "~{default_location}/~{outputname}-W~{window}-G~{gap_size}-islands-summary"
        File? fdrisland = "~{default_location}/~{outputname}-W~{window}-G~{gap_size}-FDR0.01-island.bed"
    }
}
