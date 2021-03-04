version 1.0

task sortbed {
    input {
        File bedfile

        String outputfile = basename(bedfile,'\.bed')+ '.sorted.bed'

        Int memory_gb = 10
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<
        sort-bed \
            ~{bedfile} \
            > ~{outputfile}
    >>> 
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'quay.io/biocontainers/bedops:2.4.37--hc9558a2_0'
        cpu: ncpu
    }
    output {
        File sortbed_out = "~{outputfile}"
    }
}
