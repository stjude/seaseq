version 1.0

import "https://raw.githubusercontent.com/adthrasher/seaseq/refactor/workflows/tasks/util.wdl"

workflow visualization {
    input {
        File? xlsfile
        File wigfile
        File chromsizes
        Boolean control = false
        String default_location = "Coverage_files"
    }
    
    if ( defined(xlsfile) ) {
        String string_xlsfile = "" #buffer to allow optionality
        File xlsfile_ = select_first([xlsfile, string_xlsfile])
        call util.normalize {
            input:
                wigfile=wigfile,
                control=control,
                xlsfile=xlsfile_,
                default_location=default_location
        }
    }

    File processed_wigfile = select_first([normalize.norm_wig, wigfile])

    call wigtobigwig {
        input:
            chromsizes=chromsizes,
            wigfile=processed_wigfile,
            default_location=default_location
    }
    call igvtdf {
        input:
            wigfile=processed_wigfile,
            chromsizes=chromsizes,
            default_location=default_location
    }
    
    output {
        File bigwig = wigtobigwig.bigwig
        File? norm_wig = normalize.norm_wig
        File tdffile = igvtdf.tdffile
    }
    
}
task wigtobigwig {
    input {
        File wigfile
        File chromsizes
        String default_location = "Coverage_files"

        String outputfile = sub(basename(wigfile),'\.wig\.gz', '.bw')

        Int memory_gb = 5
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<
        mkdir -p ~{default_location} && cd ~{default_location}

        wigToBigWig \
            -clip \
            ~{wigfile} \
            ~{chromsizes} \
            ~{outputfile}
    >>> 
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'abralab/wigtobigwig:v4'
        cpu: ncpu
    }
    output {
        File bigwig = "~{default_location}/~{outputfile}"
    }
}

task igvtdf {
    input {
        File wigfile
        File chromsizes
        String default_location = "Coverage_files"

        String outputfile = sub(basename(wigfile),'\.wig\.gz', '.tdf')

        Int memory_gb = 5
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<
        mkdir -p ~{default_location} && cd ~{default_location}
        ln -s ~{chromsizes} genome.chrom.sizes

        igvtools \
            toTDF \
            ~{wigfile} \
            ~{outputfile} \
            genome.chrom.sizes
    >>> 
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'abralab/igvtools:v2.8.2'
        cpu: ncpu
    }
    output {
        File tdffile = "~{default_location}/~{outputfile}"
    }
}
