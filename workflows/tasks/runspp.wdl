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
        run_spp.R \
            -c=~{basename(bamfile)} \
            ~{true="-savp" false="" crosscorr} \
            -out=~{outputfile}
    >>> 
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'madetunj/spp:v1.16.0'
        cpu: ncpu
    }
    output {
        File spp_out = "~{outputfile}"
    }
}
