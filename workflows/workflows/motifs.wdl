version 1.0

import "https://raw.githubusercontent.com/stjude/seaseq/master/workflows/tasks/bedtools.wdl"
#import "/home/madetunj/EDITS/seaseq/workflows/tasks/bedtools.wdl"

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

    call process_motif_folder {
        input :
            inputfile=bedfasta.fastafile
    }

    call ame {
        input :
            motif_databases=motif_databases,
            fastafile=bedfasta.fastafile,
            folder_output = select_first(process_motif_folder.placeholder_output),
            default_location=default_location + '/' + basename(select_first(process_motif_folder.placeholder_output)) + '-ame_out',
    }

    call meme {
        input :
            fastafile=bedfasta.fastafile,
            default_location=default_location,
            folder_output = select_first(process_motif_folder.placeholder_output)
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
        File folder_output

        Int memory_gb = 5
        Int max_retries = 1
        Int ncpu = 1

        String outputfolder = basename(folder_output) + '-meme_out'
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
        docker: 'abralab/memesuite:v5.3.3'
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
        File folder_output

        Int memory_gb = 10
        Int max_retries = 1
        Int ncpu = 1

        String default_location = "MOTIF_files" + '/' + basename(folder_output) + '-ame_out'
    }
    command <<<
        mkdir -p ~{default_location} && cd ~{default_location}

        ame \
            -oc ./ \
            ~{fastafile} \
            ~{sep=' ' motif_databases}

       gzip sequences.tsv

    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'abralab/memesuite:v5.3.3'
        cpu: ncpu
    }
    output {
        File? ame_tsv = "~{default_location}/ame.tsv"
        File? ame_html = "~{default_location}/ame.html"
        File? ame_seq = "~{default_location}/sequences.tsv.gz"
    }

}

task process_motif_folder {

    input {
        File inputfile
        String parser = basename(inputfile,'.fa')
    }

    command <<<
        old_name=~{parser}
        output_name=${old_name##*.}
        mkdir -p $old_name
        touch ~{parser}/$output_name
    >>>

    output {
        Array[File?] placeholder_output = glob("~{parser}/*")
    }

}
