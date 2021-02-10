version 1.0

task macs {

    input {
        File bamfile

        Int memory_gb = 10
        Int max_retries = 1
        Int ncpu = 1

        String pvalue = '1e-9'
        String keep_dup = 'auto'
        Int space = 50
        Int shiftsize = 200
        Boolean wiggle = true
        Boolean single_profile = true
        Boolean nomodel = false

        String prefix_location = "PEAKS_files/NARROW_peaks"
        String default_location = if (nomodel) then prefix_location + '/' + basename(bamfile,'\.bam')  + '-nm' else prefix_location + '/' + basename(bamfile,'\.bam')  + '-p' + sub(pvalue,'^.*-','') + '_kd-' + keep_dup
        String output_name = if (nomodel) then basename(bamfile,'\.bam') + '_nm' else basename(bamfile,'\.bam') + '_p' + sub(pvalue,'^.*-','') + '_kd-' + keep_dup

    }
    command <<<

        mkdir -p ~{default_location}

        if [ "~{nomodel}" == 'true' ]; then
            macs14 \
                -t ~{bamfile} \
                --space=~{space} \
                --shiftsize=~{shiftsize} \
                ~{true="-w" false="" wiggle} \
                ~{true="-S" false="" single_profile} \
                ~{true="--nomodel" false="" nomodel} \
                -n ~{output_name} \
                && mv \
                ~{output_name}* \
                ~{default_location}
        else
            macs14 \
                -t ~{bamfile} \
                -p ~{pvalue} \
                --keep-dup=~{keep_dup} \
                --space=~{space} \
                ~{true="-w" false="" wiggle} \
                ~{true="-S" false="" single_profile} \
                -n ~{output_name} \
                && mv \
                ~{output_name}* \
                ~{default_location}
        fi
    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'madetunj/macs:v1.4.2'
        cpu: ncpu
    }
    output {
        File peakbedfile = "~{default_location}\/~{output_name}_peaks.bed"
        File peakxlsfile = "~{default_location}\/~{output_name}_peaks.xls"
        File summitsfile = "~{default_location}\/~{output_name}_summits.bed"
        File wigfile = "~{default_location}\/~{output_name}_MACS_wiggle\/treat\/~{output_name}_treat_afterfiting_all.wig.gz"
    }
}
