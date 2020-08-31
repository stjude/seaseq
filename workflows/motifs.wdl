version 1.0

import "https://raw.githubusercontent.com/stjude/seaseq/wdl-workflows/tasks/bedtools.wdl"

workflow motifs {
    input {
        File bedfile
        File reference 
        File reference_index
        Array[File]+ motif_databases
        String default_location = "MOTIF_files"
    }

    call bedtools.bedfasta {
        input :
            reference=reference,
            reference_index=reference_index,
            bedfile=bedfile
    }

    call ame {
        input :
            motif_databases=motif_databases,
            fastafile=bedfasta.fastafile,
            default_location=default_location
    }

    call meme {
        input :
            fastafile=bedfasta.fastafile,
            default_location=default_location
    }

    output {
        File? ame_tsv = ame.ame_tsv
        File? ame_html = ame.ame_html
        File? ame_seq = ame.ame_seq
        File meme_out = meme.outputdir 
        File meme_summary = meme.meme_summary
    }
}

task meme {

    input {
        File fastafile
        Boolean spamo_skip = false
        Boolean fimo_skip = false
        String default_location = "MOTIF_files"

        Int memory_gb = 5
        Int max_retries = 1
        Int ncpu = 1

        String outputfolder = "bklist" + sub(basename(fastafile,'.fa'),'^.*bklist','') + '-meme_out'
    }
    command <<<
        mkdir -p ~{default_location} && cd ~{default_location}

        meme-chip \
            ~{true="-spamo-skip" false="" spamo_skip} \
            ~{true="-fimo-skip" false="" fimo_skip} \
            -oc ~{outputfolder} \
            ~{fastafile}

       zip -9r ~{outputfolder}.zip ~{outputfolder}
       cp ~{outputfolder}/summary.tsv ~{outputfolder}-summary.tsv
    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'madetunj/memesuite:v5.1.1'
        cpu: ncpu
    }
    output {
        File outputdir = "~{default_location}/~{outputfolder}.zip"
        File meme_summary = "~{default_location}/~{outputfolder}-summary.tsv"
    }
}

task ame {

    input {
        File fastafile
        Array[File]+ motif_databases
        String default_location = "MOTIF_files"

        Int memory_gb = 10
        Int max_retries = 1
        Int ncpu = 1

        String outputfolder = "bklist" + sub(basename(fastafile,'.fa'),'^.*bklist','') + '-ame_out'
    }
    command <<<
        mkdir -p ~{default_location} && cd ~{default_location}

        ame \
            -oc ~{outputfolder} \
            ~{fastafile} \
            ~{sep=' ' motif_databases}
       gzip ~{outputfolder}/sequences.tsv
    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'madetunj/memesuite:v5.1.1'
        cpu: ncpu
    }
    output {
        File? ame_tsv = "~{default_location}/~{outputfolder}/ame.tsv"
        File? ame_html = "~{default_location}/~{outputfolder}/ame.html"
        File? ame_seq = "~{default_location}/~{outputfolder}/sequences.tsv.gz"
    }
}
