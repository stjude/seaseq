version 1.0
# Bedtools tools

task intersect {
    input {
        File fileA
        File fileB
        Boolean nooverlap = false
        Boolean countoverlap = false
        Boolean sorted = false

        String outputfile = sub(basename(fileA), '\.b..$', '')
        String suffixname = if (nooverlap) then '.bklist.bam' else '.sorted.bed'
        String default_location = "."
        
        Int memory_gb = 10
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<
        mkdir -p ~{default_location} && cd ~{default_location}

        intersectBed \
            ~{true="-v" false="" nooverlap} \
            -a ~{fileA} \
            -b ~{fileB} \
            ~{true="-c" false="" countoverlap} \
            ~{true="-sorted" false="" sorted} \
            > ~{outputfile}~{suffixname}
    >>> 
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'madetunj/bedtools:v2.25.0'
        cpu: ncpu
    }
    output {
        File intersect_out = "~{default_location}/~{outputfile}~{suffixname}" 
    }
}

task bamtobed {
    input {
        File bamfile
        String outputfile = basename(bamfile) + "2bed.bed"

        Int memory_gb = 5
        Int max_retries = 1
        Int ncpu = 1
    }
    command {
        bamToBed \
            -i ~{bamfile} \
            > ~{outputfile}
    }
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'madetunj/bedtools:v2.25.0'
        cpu: ncpu
    }
    output {
        File bedfile = "~{outputfile}"
    }
}

task bedfasta {
    input {
        File bedfile
        String outputfile = basename(bedfile,'.bed') + ".fa"
        File reference
        File reference_index

        Int memory_gb = 5
        Int max_retries = 1
        Int ncpu = 1
    }
    command {
        ln -s ~{reference} ~{basename(reference)}
        ln -s ~{reference_index} ~{basename(reference_index)}
        bedtools \
            getfasta \
            -fi ~{basename(reference)} \
            -bed ~{bedfile} \
            -fo ~{outputfile}
    }
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'madetunj/bedtools:v2.25.0'
        cpu: ncpu
    }
    output {
        File fastafile = "~{outputfile}"
    }
}
