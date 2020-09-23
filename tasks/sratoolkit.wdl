version 1.0

task fastqdump {
    input {
        String sraid
        Int memory_gb = 5
        Int max_retries = 1
        Int ncpu = 1
    }
    command {
        fastq-dump --gzip ~{sraid}
    }
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'madetunj/sratoolkit:v2.9.6'
        cpu: ncpu
    }
    output {
        File srafastq = "~{sraid}.fastq.gz"
    }
}
