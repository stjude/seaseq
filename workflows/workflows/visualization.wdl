version 1.0

import "https://raw.githubusercontent.com/stjude/seaseq/master/workflows/tasks/util.wdl"

workflow visualization {
    input {
        File xlsfile
        File wigfile
        File chromsizes
        String default_location = "PEAKDisplay_files"
    }
    
    call util.normalize {
        input:
            wigfile=wigfile,
            xlsfile=xlsfile,
            default_location=default_location
    }
    call wigtobigwig {
        input:
            chromsizes=chromsizes,
            wigfile=normalize.norm_wig,
            default_location=default_location
    }
    call igvtdf {
        input:
            wigfile=normalize.norm_wig,
            default_location=default_location
    }
    
    output {
        File bigwig = wigtobigwig.bigwig
        File norm_wig = normalize.norm_wig
        File tdffile = igvtdf.tdffile
    }
    
}
task wigtobigwig {
    input {
        File wigfile
        File chromsizes
        String default_location = "PEAKDisplay_files"

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
        docker: 'madetunj/wigtobigwig:v4'
        cpu: ncpu
    }
    output {
        File bigwig = "~{default_location}/~{outputfile}"
    }
}

task igvtdf {
    input {
        File wigfile
        String default_location = "PEAKDisplay_files"

        String genome = "hg19"

        String outputfile = sub(basename(wigfile),'\.wig\.gz', '.tdf')

        Int memory_gb = 5
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<
        mkdir -p ~{default_location} && cd ~{default_location}

        igvtools \
            toTDF \
            ~{wigfile} \
            ~{outputfile} \
            ~{genome}
    >>> 
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'madetunj/igvtools:v2.8.2'
        cpu: ncpu
    }
    output {
        File tdffile = "~{default_location}/~{outputfile}"
    }
}
