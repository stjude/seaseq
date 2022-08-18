version 1.0
# SAMtools

task indexstats {
    input {
        File bamfile
        String outputfile = basename(bamfile) + ".bai"
        String flagstat = sub(basename(bamfile),".bam$", "-flagstat.txt")
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
        docker: 'ghcr.io/stjude/abralab/samtools:v1.9'
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
        String outputfile = sub(basename(bamfile),".bam$", ".rmdup.bam")
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
        docker: 'ghcr.io/stjude/abralab/samtools:v1.9'
        cpu: ncpu
    }
    output {
        File mkdupbam = "~{default_location}/~{outputfile}"
    }
}

task viewsort {
    input {
        File samfile
        String outputfile = basename(sub(samfile,'.sam$','.sorted.bam'))
        String fixmatefile = basename(sub(samfile,'.sam$','.fixmate.bam'))
        String default_location = "BAM_files"
        Boolean paired_end = false

        Int memory_gb = 5
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<
        mkdir -p ~{default_location} && cd ~{default_location}

        if [ "~{paired_end}" == 'true' ]; then
            awk -F\\t 'BEGIN{j=0}{j++}{if(NF>5 && j%2==0){ \
                printf "%s_%.0f\t", $1, j-1 } else if(NF>5 && j%2==1){ \
                printf "%s_%.0f\t", $1, j } else { printf $1 "\t";j=0 } \
                for(i=2;i<=NF;i++){ printf "%s\t", $i}; printf "\n" }' \
                ~{samfile} > ~{basename(sub(samfile,'.sam','.renamed.sam'))}

            samtools fixmate -m \
                ~{basename(sub(samfile,'.sam','.renamed.sam'))} \
                ~{fixmatefile}
            samtools sort \
                ~{fixmatefile} \
                -o ~{outputfile}
        else
            samtools sort \
                ~{samfile} \
                -o ~{outputfile}
        fi

    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'ghcr.io/stjude/abralab/samtools:v1.9'
        cpu: ncpu
    }
    output {
        File sortedbam = "~{default_location}/~{outputfile}"
        File? fixmatebam = "~{default_location}/~{fixmatefile}"
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
        docker: 'ghcr.io/stjude/abralab/samtools:v1.9'
        cpu: ncpu
    }
    output {
        File faidx_file = "~{sub(basename(reference),'.gz','')}.fai"
        File chromsizes = "~{sub(basename(reference),'.gz','')}.tab"
    }
}

task mergebam {
    input {
        Array[File] bamfiles
        Array[File] metricsfiles
        String outputfile = 'AllMapped.' + length(bamfiles) + '_merge.bam'
        String fixmatefile = 'AllMapped.' + length(bamfiles) + '_merge.fixmate.bam'
        String default_location = "BAM_files"
        Boolean paired_end = false

        Int memory_gb = 5
        Int max_retries = 1
        Int ncpu = 20
    }
    command <<<
        mkdir -p ~{default_location} && cd ~{default_location}

        total_readlength=0
        for each in $(ls -1 ~{sep=' ' metricsfiles}); do
            readlength=$(tail -n 1 $each | awk '{print $4}');
            total_readlength=$(echo "$readlength + $total_readlength" | bc)
        done
        echo "$total_readlength/~{length(metricsfiles)}" | bc > average_readlength.txt

        cat average_readlength.txt

        if [ "~{paired_end}" == 'true' ]; then
            samtools merge \
                --threads ~{ncpu} \
                ~{fixmatefile} \
                ~{sep=' ' bamfiles}
            samtools sort \
                ~{fixmatefile} \
                -o ~{outputfile}
        else
            samtools merge \
                --threads ~{ncpu} \
                ~{outputfile} \
                ~{sep=' ' bamfiles}
        fi
    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'ghcr.io/stjude/abralab/samtools:v1.9'
        cpu: ncpu
    }
    output {
        File mergebam = "~{default_location}/~{outputfile}"
        Int avg_readlength = read_int("~{default_location}/average_readlength.txt")
        File? fixmatemergebam = "~{default_location}/~{fixmatefile}"
    }
}
