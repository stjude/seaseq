version 1.0

task sicer {

    input {
        File bedfile
        String default_location = "PEAKS_files/BROAD_peaks"

        String species = "hg19"
        Int redundancy = 1
        Int window = 200
        Int fragment_size = 150
        Float genome_fraction = 0.86
        Int gap_size = 200
        Int evalue = 100

        String outputname = basename(bedfile,'.bed')

        Int memory_gb = 15
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<
        mkdir -p ~{default_location}

        sicer \
            -t ~{bedfile} \
            -s ~{species} \
            -rt ~{redundancy} \
            -w ~{window} \
            -f ~{fragment_size} \
            -egf ~{genome_fraction} \
            -g ~{gap_size} \
            -e ~{evalue}

        gzip *wig
        mv *W200* ~{default_location}
    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'madetunj/sicer:v1.0.2'
        cpu: ncpu
    }
    output {
        File scoreisland = "~{default_location}/~{outputname}-W200-G200.scoreisland"
        File wigfile = "~{default_location}/~{outputname}-W200-normalized.wig.gz"
    }
}
