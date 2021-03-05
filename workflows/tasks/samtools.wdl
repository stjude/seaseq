version 1.0
# SAMtools

task indexstats {
    input {
        File bamfile
        String outputfile = basename(bamfile) + ".bai"
        String flagstat = sub(basename(bamfile),"\.bam$", "-flagstat.txt")
        String default_location = "BAM_files"

        Int memory_gb = 5
        Int max_retries = 1
        Int ncpu = 1
    }
    command {
        mkdir -p ~{default_location} && cd ~{default_location}

        ln -s ~{bamfile} ~{basename(bamfile)}

        samtools flagstat ~{bamfile} > ~{flagstat}

        samtools index ~{basename(bamfile)}
    }
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'quay.io/biocontainers/samtools:1.9--h10a08f8_12'
        cpu: ncpu
    }
    output {
        File indexbam = "~{default_location}/~{outputfile}"
        File flagstats = "~{default_location}/~{flagstat}"
    }
}

task markdup {
    input {
        File bamfile
        String outputfile = sub(basename(bamfile),"\.bam$", ".rmdup.bam")
        String default_location = "BAM_files"

        Int memory_gb = 5
        Int max_retries = 1
        Int ncpu = 1
    }
    command {
        mkdir -p ~{default_location} && cd ~{default_location}

        samtools markdup \
            -r -s \
            ~{bamfile} \
            ~{outputfile}
    }
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'quay.io/biocontainers/samtools:1.9--h10a08f8_12'
        cpu: ncpu
    }
    output {
        File mkdupbam = "~{default_location}/~{outputfile}"
    }
}

task viewsort {
    input {
        File samfile
        String outputfile = basename(sub(samfile,'\.sam$','\.sorted.bam'))
        String default_location = "BAM_files"

        Int memory_gb = 5
        Int max_retries = 1
        Int ncpu = 1
    }
    command {
        mkdir -p ~{default_location} && cd ~{default_location}

        samtools view -b \
            ~{samfile} \
            > ~{sub(samfile,'\.sam$','\.bam')}

        samtools sort \
           ~{sub(samfile,'\.sam$','\.bam')} \
           -o ~{outputfile}
    }
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'quay.io/biocontainers/samtools:1.9--h10a08f8_12'
        cpu: ncpu
    }
    output {
        File sortedbam = "~{default_location}/~{outputfile}"
    }
}

task faidx {
    input {
	File reference

        Int memory_gb = 5
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<
        ln -s ~{reference} ~{basename(reference)}
        samtools faidx ~{basename(reference)} -o ~{basename(reference)}.fai
        cut -f1,2 ~{basename(reference)}.fai > ~{basename(reference)}.tab
    >>>
    runtime {
        memory: memory_gb + " GB"
        maxRetries: max_retries
        docker: 'quay.io/biocontainers/samtools:1.9--h10a08f8_12'
        cpu: ncpu
    }
    output {
	File faidx_file = "~{basename(reference)}.fai"
        File chromsizes = "~{basename(reference)}.tab"
    }
}
