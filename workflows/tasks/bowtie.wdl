version 1.0

task bowtie {
    input {
        File fastqfile
        File? metricsfile
        Array[File]+ index_files

        String outputfile = sub(basename(fastqfile),'.fastq.gz|.fq.gz','.sam')
        
        Int? read_length = 75
        Int limit_alignments = 2
        Int good_alignments = 2
        Boolean best_alignments = true

        Int additional_memory_gb = 10
        Int max_retries = 1
        Int ncpu = 20
    }

    Int memory_gb = ceil(size(fastqfile)) + ceil(size(index_files)) + additional_memory_gb

    command <<<
        if [ -f "~{metricsfile}" ]; then
            readlength=$(tail -n 1 ~{metricsfile} | awk '{print $4}');
            echo "Metrics file with readlength " $readlength
        else
            readlength=~{read_length}
        fi

        bowtie \
            -l $readlength \
            -p ~{ncpu} \
            -k ~{good_alignments} \
            -m ~{limit_alignments} \
            ~{true="--best" false="" best_alignments} \
            -S \
            ~{sub(index_files[0], "(.rev)?.[0-9].ebwt$", "")} \
            ~{fastqfile} \
            > ~{outputfile}
    >>>
    runtime {
        memory: memory_gb + " GB"
        maxRetries: max_retries
        docker: 'abralab/bowtie:v1.2.3'
        cpu: ncpu
    }
    output {
        File samfile = "~{outputfile}"
        #Array[File?] samfile = glob("*.sam")
    }
}

task index {
    input {
        File reference

        Int memory_gb = 20
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<
        if [[ "~{reference}" == *"gz" ]]; then
            gunzip -c ~{reference} > ~{sub(basename(reference),'.gz','')}
        else
           ln -s ~{reference} ~{sub(basename(reference),'.gz','')}
        fi

        bowtie-build --threads ~{ncpu} ~{sub(basename(reference),'.gz','')} ~{sub(basename(reference),'.gz','')}-index
    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'abralab/bowtie:v1.2.3'
        cpu: ncpu
    }
    output {
        Array[File] bowtie_indexes = glob("~{sub(basename(reference),'.gz','')}-index*")
    }
}
