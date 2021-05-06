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
        docker: 'abralab/samtools:v1.9'
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
        docker: 'abralab/samtools:v1.9'
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
        docker: 'abralab/samtools:v1.9'
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
        if [[ "~{reference}" == *"gz" ]]; then
            gunzip -c ~{reference} > ~{sub(basename(reference),'.gz','')}
        else
           ln -s ~{reference} ~{sub(basename(reference),'.gz','')}
        fi

        samtools faidx ~{sub(basename(reference),'.gz','')} -o ~{sub(basename(reference),'.gz','')}.fai
        cut -f1,2 ~{sub(basename(reference),'.gz','')}.fai > ~{sub(basename(reference),'.gz','')}.tab
    >>>
    runtime {
        memory: memory_gb + " GB"
        maxRetries: max_retries
        docker: 'abralab/samtools:v1.9'
        cpu: ncpu
    }
    output {
	File faidx_file = "~{sub(basename(reference),'.gz','')}.fai"
        File chromsizes = "~{sub(basename(reference),'.gz','')}.tab"
    }
}
