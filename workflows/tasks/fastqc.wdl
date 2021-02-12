version 1.0

task fastqc {
    input {
        File inputfile
        String prefix = sub(basename(inputfile),'\.f.*q\.gz','')
        String default_location="QC_files/FASTQC"

        Int memory_gb = 5
        Int max_retries = 1
        Int ncpu = 1
    }
    command {
        ln -s ~{inputfile} ~{sub(basename(inputfile),'\.bam$','.bam.bam')}

        mkdir -p ~{default_location}

        fastqc \
            -o ~{default_location} \
            ~{sub(basename(inputfile),'\.bam$','.bam.bam')}
    }
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'madetunj/fastqc:v0.11.9'
        cpu: ncpu
    }
    output {
        File htmlfile = "~{default_location}/~{prefix}_fastqc.html"
        File zipfile = "~{default_location}/~{prefix}_fastqc.zip"
    }
}
