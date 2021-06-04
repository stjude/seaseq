version 1.0

task bowtie2 {
    input {
        File fastqfile
        File? metricsfile
        Array[File]+ index_files

        String outputfile = sub(basename(fastqfile),'\.f.*q\.gz','.sam')
        
        Int? read_length = 75
        Int limit_alignments = 2
        Int good_alignments = 2
        Boolean best_alignments = true

        Int memory_gb = 10
        Int max_retries = 1
        Int ncpu = 20
    }
    command <<<
        if [ -f "~{metricsfile}" ]; then
            readlength=$(tail -n 1 ~{metricsfile} | awk '{print $4}');
        else
            readlength=~{read_length}
        fi

        bowtie2 \
            -p ~{ncpu} \
            -x ~{sub(index_files[0], "(\.rev)?\.[0-9]\.bt2$", "")} \
            -U ~{fastqfile} \
            -S ~{outputfile}
    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'abralab/bowtie2:v2.4.4'
        cpu: ncpu
    }
    output {
        Array[File?] samfile = glob("*.sam")
    }
}

task index {
    input {
        File reference

        Int memory_gb = 10
        Int max_retries = 1
        Int ncpu = 20
    }
    command <<<
        if [[ "~{reference}" == *"gz" ]]; then
            gunzip -c ~{reference} > ~{sub(basename(reference),'.gz','')}
        else
           ln -s ~{reference} ~{sub(basename(reference),'.gz','')}
        fi

        bowtie2-build --threads ~{ncpu} ~{sub(basename(reference),'.gz','')} ~{sub(basename(reference),'.gz','')}-index
    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'abralab/bowtie2:v2.4.4'
        cpu: ncpu
    }
    output {
        Array[File] bowtie_indexes = glob("~{sub(basename(reference),'.gz','')}-index*")
    }
}
