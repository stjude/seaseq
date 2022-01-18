version 1.0

task runspp {
    input {
        File bamfile
        File? control
        Boolean crosscorr = true

        String outputfile = basename(bamfile,'.bam')+ '-spp.out'

        Int memory_gb = 10
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<
        ln -s ~{bamfile} ~{basename(bamfile)}
        ln -s ~{control} control.bam

        mappedreads=$(samtools view -F 0x0204 ~{basename(bamfile)} | wc -l)
        if [ $mappedreads -gt 0 ]; then
            run_spp.R \
                -c=~{basename(bamfile)} \
                ~{if defined(control) then "-i=control.bam" else ""} \
                ~{true="-savp" false="" crosscorr} \
                -out=~{outputfile}
        else
            touch ~{outputfile}
        fi
    >>> 
    runtime {
        continueOnReturnCode: [0, 1]
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'abralab/spp:v1.16.0'
        cpu: ncpu
    }
    output {
        File? spp_out = "~{outputfile}"
    }
}
