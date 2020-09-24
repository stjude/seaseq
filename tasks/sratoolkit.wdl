version 1.0

task fastqdump {
    input {
        String? sra_id

        Int memory_gb = 5
        Int max_retries = 1
        Int ncpu = 20
    }
    command {
        pfastq-dump -t 20 --gzip -s ~{sra_id} -O ./
    }
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'madetunj/sratoolkit:v2.9.6'
        cpu: ncpu
    }
    output {
        File fastqfile = "~{sra_id}.fastq.gz"
    }
}
