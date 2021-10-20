version 1.0

task fastqdump {
    input {
        String? sra_id
        Boolean cloud=false
        Boolean paired=false
        Int memory_gb = 5
        Int max_retries = 1
        Int ncpu = 20
    }
    command {
        threads=~{ncpu}
        if [ "~{cloud}" == 'true' ]; then
            output=$(sra-stat --meta --quick ~{sra_id})
            total=0
            for line in $output; do
                value=$(echo $line | cut -d '|' -f 3 | cut -d ':' -f 1)
                total=$(($total + $value))
            done
            if [ "$total" -ge 50000000 ]; then
                if [ "$threads" -ge 10 ]; then threads=10; fi
                echo "number of threads changed to 10 for dnanexus compatibility"
            fi
        fi
        pfastq-dump \
            -t $threads \
            --gzip \
            ~{true="--split-files" false="" paired} \
            -s ~{sra_id} -O ./
    }
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'abralab/sratoolkit:v2.9.6'
        cpu: ncpu
    }
    output {
        Array[File] fastqfile = glob("~{sra_id}.fastq.gz")
    }
}
