version 1.0

task runspp {
    input {
        File bamfile
        Boolean crosscorr = true

        String outputfile = basename(bamfile,'\.bam')+ '-spp.out'

        Int memory_gb = 10
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<
        ln -s ~{bamfile} ~{basename(bamfile)}

        mappedreads=$(samtools view -F 0x0204 ~{basename(bamfile)} | wc -l)
        if [ $mappedreads -gt 0 ]; then
            run_spp.R \
                -c=~{basename(bamfile)} \
                ~{true="-savp" false="" crosscorr} \
                -out=~{outputfile}
        else
            touch ~{outputfile}
        fi
    >>> 
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'abralab/spp:v1.16.0'
        cpu: ncpu
    }
    output {
        File spp_out = "~{outputfile}"
    }
}
