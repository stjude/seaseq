version 1.0

task macs {

    input {
        File bamfile
        File? control

        Int memory_gb = 25
        Int max_retries = 1
        Int ncpu = 1

        String pvalue = '1e-9'
        String keep_dup = 'auto'
        String egs = 'hs'
        Int space = 50
        Int shiftsize = 200
        Boolean wiggle = true
        Boolean single_profile = true
        Boolean nomodel = false

        String prefix_location = "PEAKS_files/NARROW_peaks"
        String coverage_location = "COVERAGE_files/NARROW_peaks"
        String default_location = if (nomodel) then prefix_location + '/' + basename(bamfile,'.bam')  + '-nm' else prefix_location + '/' + basename(bamfile,'.bam')  + '-p' + sub(pvalue,'^.*-','') + '_kd-' + keep_dup
        String output_name = if (nomodel) then basename(bamfile,'.bam') + '_nm' else basename(bamfile,'.bam') + '_p' + sub(pvalue,'^.*-','') + '_kd-' + keep_dup

    }
    command <<<

        mkdir -p ~{default_location} ~{coverage_location}/control

        if [ "~{nomodel}" == 'true' ]; then
            macs14 \
                -t ~{bamfile} \
                ~{if defined(control) then "-c" + control else ""} \
                --space=~{space} \
                --gsize=~{egs} \
                --shiftsize=~{shiftsize} \
                ~{true="-w" false="" wiggle} \
                ~{true="-S" false="" single_profile} \
                ~{true="--nomodel" false="" nomodel} \
                -n ~{output_name}
        else
            macs14 \
                -t ~{bamfile} \
                ~{if defined(control) then "-c" + control else ""} \
                -p ~{pvalue} \
                --keep-dup=~{keep_dup} \
                --space=~{space} \
                --gsize=~{egs} \
                ~{true="-w" false="" wiggle} \
                ~{true="-S" false="" single_profile} \
                -n ~{output_name}
        fi

        mv ~{output_name}_MACS_wiggle/treat/~{output_name}_treat_afterfiting_all.wig.gz ~{coverage_location}
        
        mv ~{output_name}* ~{default_location}

        #move coverage files to new location
        #mv ~{default_location}/~{output_name}_MACS_wiggle/treat/~{output_name}_treat_afterfiting_all.wig.gz ~{coverage_location}
        if [ -f "~{control}" ]; then
            mv ~{default_location}/~{output_name}_MACS_wiggle/control/~{output_name}_control_afterfiting_all.wig.gz ~{coverage_location}/control/
        fi
        
    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'ghcr.io/stjude/abralab/macs:v1.4.2'
        cpu: ncpu
    }
    output {
        File peakbedfile = "~{default_location}/~{output_name}_peaks.bed"
        File peakxlsfile = "~{default_location}/~{output_name}_peaks.xls"
        File summitsfile = "~{default_location}/~{output_name}_summits.bed"
        File? negativepeaks = "~{default_location}/~{output_name}_negative_peaks.xls"
        File wigfile = "~{coverage_location}/~{output_name}_treat_afterfiting_all.wig.gz"
        File? ctrlwigfile = "~{coverage_location}/control/~{output_name}_control_afterfiting_all.wig.gz"
    }
}
