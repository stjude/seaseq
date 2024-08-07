version 1.0

task bowtie {
    input {
        File fastqfile
        File? fastqfile_R2
        File? metricsfile
        Array[File]+ index_files
        String? prefix
        String outputfile = if (defined(prefix)) then select_first([prefix,basename(fastqfile)]) + '.sam' else if (defined(fastqfile_R2)) then sub(basename(fastqfile),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','.sam') else sub(basename(fastqfile),'.fastq.gz|.fq.gz','.sam')
        Int read_length = 75
        Int insert_size = 600
        Int limit_alignments = 2
        Int good_alignments = 2
        Boolean best_alignments = true
        String strandedness = 'fr'
        String stranded_m = if strandedness=='fr' then '--fr'
                        else if strandedness=='rf' then '--rf'
                        else if strandedness=='ff' then '--ff'
                        else '--fr'

        Int additional_memory_gb = 10
        Int max_retries = 1
        Int ncpu = 20
    }

    Int memory_gb = ceil(size(fastqfile, "GiB")) + ceil(size(index_files, "GiB")) + additional_memory_gb

    command <<<
        if [ -f "~{metricsfile}" ]; then
            readlength=$(tail -n 1 ~{metricsfile} | awk '{print $6}');
            if [ $readlength -lt 28 ]; then readlength=28; fi
            echo "Metrics file with readlength " $readlength
        else
            readlength=~{read_length}
        fi
        
        echo $readlength
        if [ -f "~{fastqfile_R2}" ]; then
            bowtie \
                --chunkmbs=256 \
                -l $readlength \
                -p ~{ncpu} \
                -k ~{good_alignments} \
                -m ~{limit_alignments} \
                -X ~{insert_size} \
                ~{stranded_m} \
                ~{true="--best" false="" best_alignments} \
                -S \
                ~{sub(index_files[0], "(.rev)?.[0-9].ebwt$", "")} \
                -1 ~{fastqfile} \
                -2 ~{fastqfile_R2} \
                > ~{outputfile}

        else
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
                
        fi
    >>>
    runtime {
        memory: memory_gb + " GB"
        maxRetries: max_retries
        docker: 'ghcr.io/stjude/abralab/bowtie:v1.2.3'
        cpu: ncpu
    }
    output {
        File samfile = "~{outputfile}"
        #Array[File?] samfile = glob("*.sam")
    }
}
task spikein_SE {
    input {
        File fastqfile
        File? metricsfile
        Array[File]+ index_files
        String? prefix
        String outputfile = if (defined(prefix)) then select_first([prefix,basename(fastqfile)]) + '.unaligned' else sub(basename(fastqfile),'.fastq.gz|.fq.gz','.unaligned')
        String nameoffile = basename(fastqfile,'.gz')
        Int read_length = 75
        Int insert_size = 600
        Boolean best_alignments = true
        String strandedness = 'fr'
        String stranded_m = if strandedness=='fr' then '--fr'
                        else if strandedness=='rf' then '--rf'
                        else if strandedness=='ff' then '--ff'
                        else '--fr'
        String default_location = "."
        Int additional_memory_gb = 10
        Int max_retries = 1
        Int ncpu = 20
    }

    Int memory_gb = ceil(size(fastqfile, "GiB")) + ceil(size(index_files, "GiB")) + additional_memory_gb

    command <<<
        mkdir -p ~{default_location}

        if [ -f "~{metricsfile}" ]; then
            readlength=$(tail -n 1 ~{metricsfile} | awk '{print $6}');
            if [ $readlength -lt 28 ]; then readlength=28; fi
            echo "Metrics file with readlength " $readlength
        else
            readlength=~{read_length}
        fi
        
        echo $readlength
        bowtie \
            -l $readlength \
            -p ~{ncpu} \
            --un ~{nameoffile} \
            ~{true="--best" false="" best_alignments} \
            -S \
            ~{sub(index_files[0], "(.rev)?.[0-9].ebwt$", "")} \
            ~{fastqfile} \
            1>output.sam 2>~{default_location}/~{outputfile}.out
        
        gzip ~{nameoffile}
    >>>
    runtime {
        memory: memory_gb + " GB"
        maxRetries: max_retries
        docker: 'ghcr.io/stjude/abralab/bowtie:v1.2.3'
        cpu: ncpu
    }
    output {
        File unaligned = "~{nameoffile}.gz"
        File mapping_output = "~{default_location}/~{outputfile}.out"
    }
}

task spikein_PE {
    input {
        File fastqfile
        File fastqfile_R2
        File? metricsfile
        Array[File]+ index_files
        String? prefix
        String outputfile = if (defined(prefix)) then select_first([prefix,basename(fastqfile)]) + '.unaligned' else if (defined(fastqfile_R2)) then sub(basename(fastqfile),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','.unaligned') else sub(basename(fastqfile),'.fastq.gz|.fq.gz','.unaligned')
        String nameoffile_R1 = basename(fastqfile,'.gz')
        String nameoffile_R2 = basename(fastqfile_R2,'.gz')
        Int read_length = 75
        Int insert_size = 600
        Boolean best_alignments = true
        String strandedness = 'fr'
        String stranded_m = if strandedness=='fr' then '--fr'
                        else if strandedness=='rf' then '--rf'
                        else if strandedness=='ff' then '--ff'
                        else '--fr'
        String default_location = "."

        Int additional_memory_gb = 10
        Int max_retries = 1
        Int ncpu = 20
    }

    Int memory_gb = ceil(size(fastqfile, "GiB")) + ceil(size(index_files, "GiB")) + additional_memory_gb

    command <<<
        if [ -f "~{metricsfile}" ]; then
            readlength=$(tail -n 1 ~{metricsfile} | awk '{print $6}');
            if [ $readlength -lt 28 ]; then readlength=28; fi
            echo "Metrics file with readlength " $readlength
        else
            readlength=~{read_length}
        fi
        
        echo $readlength
        bowtie \
            --chunkmbs=256 \
            -l $readlength \
            --un unaligned \
            -p ~{ncpu} \
            -X ~{insert_size} \
            ~{stranded_m} \
            ~{true="--best" false="" best_alignments} \
            -S \
            ~{sub(index_files[0], "(.rev)?.[0-9].ebwt$", "")} \
            -1 ~{fastqfile} \
            -2 ~{fastqfile_R2} \
            1>output.sam 2>~{outputfile}.out
            
        mv unaligned_1 ~{nameoffile_R1}; gzip ~{nameoffile_R1}
        mv unaligned_2 ~{nameoffile_R2}; gzip ~{nameoffile_R2}
        
    >>>
    runtime {
        memory: memory_gb + " GB"
        maxRetries: max_retries
        docker: 'ghcr.io/stjude/abralab/bowtie:v1.2.3'
        cpu: ncpu
    }
    output {
        File unaligned_1 = "~{nameoffile_R1}.gz"
        File unaligned_2 = "~{nameoffile_R2}.gz"
        File mapping_output = "~{outputfile}.out"
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
        docker: 'ghcr.io/stjude/abralab/bowtie:v1.2.3'
        cpu: ncpu
    }
    output {
        Array[File] bowtie_indexes = glob("~{sub(basename(reference),'.gz','')}-index*")
    }
}
