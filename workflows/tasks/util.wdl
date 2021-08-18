version 1.0

task basicfastqstats {
    input {
        File fastqfile
        String outputfile = sub(basename(fastqfile),'\.f.*q\.gz', '-fastq.readlength_dist.txt')
        String default_location = "QC_files/STATS"

        Int max_retries = 1
        Int ncpu = 1
    }

    Int memory_gb = ceil((size(fastqfile, "GiB") * 2) + 10)

    command <<<

        mkdir -p ~{default_location} && cd ~{default_location}

        zcat ~{fastqfile} | awk 'NR%4==2' | awk '{print length}' | sort -n > values.dat
        
        stddev=$(awk '{x+=$0;y+=$0^2}END{print sqrt(y/NR-(x/NR)^2)}' values.dat)

        median=$(awk '{ a[i++]=$1; } END { print a[int(i/2)]; }' values.dat)
        
        average=$(awk '{ sum += $1 } END { if (NR > 0) print sum / NR }' values.dat)
        
        minimum=$(head -n 1 values.dat)
        
        maximum=$(tail -n 1 values.dat)
        
        Q1=$(awk 'BEGIN{c=0} {total[c]=$1; c++;} END{print total[int(NR*0.25 - 0.5)]}' values.dat)
        
        Q3=$(awk 'BEGIN{c=0} {total[c]=$1; c++;} END{print total[int(NR*0.75 - 0.5)]}' values.dat)
        
        IQR=$(echo "$Q3-$Q1" | bc)
        
        
        echo Min.$'\t'1st Qu.$'\t'Median$'\t'Mean$'\t'3rd Qu.$'\t'Max.$'\t'StdDev.$'\t'IQR > ~{outputfile}
        echo $minimum$'\t'$Q1$'\t'$median$'\t'$average$'\t'$Q3$'\t'$maximum$'\t'$stddev$'\t'$IQR >> ~{outputfile}
        
        rm -rf values.dat

    >>> 
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'abralab/seaseq:v2.0.0'
        cpu: ncpu
    }
    output {
        File metrics_out = "~{default_location}/~{outputfile}"
    }
}


task flankbed {
    input {
        File bedfile
        String outputfile = basename(bedfile, '.bed') + '-flank' + flank + '.bed'
        String default_location = "."

        Int flank = 50

        Int memory_gb = 5
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<

        mkdir -p ~{default_location} && cd ~{default_location}
        ln -s ~{bedfile} ~{basename(bedfile)}
        echo ~{flank} > flank.rdl
        python <<CODE

        input = open("~{basename(bedfile)}",'r')
        output = open("~{outputfile}",'w')
        flank = int(open("flank.rdl").readline().rstrip())
        Lines = input.readlines()
        for line in Lines :
            all = line.rstrip("\n").split('\t')
            start = int(all[1]) - flank
            if start <= 0 : start = 1
            end = int(all[1]) + flank
            if end <= 0 : end = 1
            output.write("%s\t%d\t%d\t%s\t%s\n" %(all[0], start, end, all[3], all[4]))
        CODE
        rm -rf flank.rdl

    >>> 
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'abralab/seaseq:v2.0.0'
        cpu: ncpu
    }
    output {
        File flankbedfile = "~{default_location}/~{outputfile}"
    }
}


task summaryreport {
    input {
        File? controlqc_txt
        File? controlqc_html
        File sampleqc_txt
        File sampleqc_html
        File overallqc_txt
        File overallqc_html

        String outputfile = sub(basename(overallqc_html), 'stats.html', 'seaseq_report.html')
        String outputtxt = sub(basename(overallqc_html), 'stats.html', 'seaseq_report.txt')
        
        Int memory_gb = 5
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<

        # Printing Sample Quality Reports
        head -n 101 /usr/local/bin/scripts/seaseq_overall.header > ~{outputfile}
        cat ~{sampleqc_html} >> ~{outputfile}
        echo -e '</table></div>\n<div class="body">' >> ~{outputfile}

        echo -e 'SEAseq Report\nSEAseq Quality Statistics and Evaluation Report\n\nSample FASTQs Quality Results' > ~{outputtxt}
        cat ~{sampleqc_txt} >> ~{outputtxt}
        echo -e '\n' >> ~{outputtxt}

        if [ -f "~{controlqc_html}" ]; then
            # Printing Control Quality Reports
            echo '<h2>Control FASTQs Quality Results</h2>' >> ~{outputfile}
            cat ~{controlqc_html} >> ~{outputfile}
            echo -e '</table></div>\n<div class="body">' >> ~{outputfile}

            echo 'Control FASTQs Quality Results' >> ~{outputtxt}
            cat ~{controlqc_txt} >> ~{outputtxt}
            echo -e '\n' >> ~{outputtxt}
        fi
        
        # Printing Overall Quality Reports
        echo '<h2>Overall Quality Evaluation and Statistics Results</h2><p>' >> ~{outputfile}
        cat ~{overallqc_html} >> ~{outputfile}
        echo -e '</table></div>\n</body>\n</html>' >> ~{outputfile}

        echo -e 'Overall Quality Evaluation and Statistics Results' >> ~{outputtxt}
        cat ~{overallqc_txt} >> ~{outputtxt}
        echo -e '\n' >> ~{outputtxt}

    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'abralab/seaseq:v2.0.0'
        cpu: ncpu
    }
    output {
        File summaryhtml = "~{outputfile}"
        File summarytxt = "~{outputtxt}"
    }
}


task evalstats {
    input {
        File? sample_bambed
        File? control_bambed
        File? sppfile
        File? sample_countsfile
        File? sample_peaksxls
        File? control_countsfile
        File? control_peaksxls
        File sample_fastqczip
        File? control_fastqczip
        File? sample_bamflag
        File? sample_rmdupflag
        File? sample_bkflag
        File? control_bamflag
        File? control_rmdupflag
        File? control_bkflag
        File? fastqmetrics
        File? enhancers
        File? superenhancers

        String fastq_type = "Sample FASTQs"
        String default_location = "QC_files/STATS"
        String outputfile = sub(basename(sample_fastqczip),'_fastqc.zip', '-stats.csv')
        String outputhtml = sub(basename(sample_fastqczip),'_fastqc.zip', '-stats.html')
        String outputtext = sub(basename(sample_fastqczip),'_fastqc.zip', '-stats.txt')

        Int memory_gb = 10
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<

        mkdir -p ~{default_location}
        cd ~{default_location}

        evaluation-statistics.pl \
            -fqc ~{sample_fastqczip} \
            ~{if defined(sample_bambed) then "-b " + sample_bambed else ""} \
            ~{if defined(control_bambed) then "-cb " + control_bambed else ""} \
            ~{if defined(sppfile) then "-s " + sppfile else ""} \
            ~{if defined(sample_countsfile) then "-c " + sample_countsfile else ""} \
            ~{if defined(sample_peaksxls) then "-px " + sample_peaksxls else ""} \
            ~{if defined(control_countsfile) then "-ccount " + control_countsfile else ""} \
            ~{if defined(control_peaksxls) then "-cpx " + control_peaksxls else ""} \
            ~{if defined(sample_bamflag) then "-bamflag " + sample_bamflag else ""} \
            ~{if defined(sample_bkflag) then "-bkflag " + sample_bkflag else ""} \
            ~{if defined(sample_rmdupflag) then "-rmdupflag " + sample_rmdupflag else ""} \
            ~{if defined(control_bamflag) then "-cbamflag " + control_bamflag else ""} \
            ~{if defined(control_bkflag) then "-cbkflag " + control_bkflag else ""} \
            ~{if defined(control_rmdupflag) then "-crmdupflag " + control_rmdupflag else ""} \
            ~{if defined(control_fastqczip) then "-cfqc " + control_fastqczip else ""} \
            ~{if defined(fastqmetrics) then "-fx " + fastqmetrics else ""} \
            ~{if defined(enhancers) then "-re " + enhancers else ""} \
            ~{if defined(superenhancers) then "-rs " + superenhancers else ""} \
            -outfile ~{outputfile}

        tail -n 101 /usr/local/bin/scripts/seaseq_overall.header > ~{outputhtml}
        cat ~{outputhtml}x >> ~{outputhtml}
        sed -i "s/SEAseq ~{fastq_type} Report/SEAseq ~{fastq_type} Report/" ~{outputhtml}


    >>> 
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'abralab/seaseq:v2.0.0'
        cpu: ncpu
    }
    output {
        File statsfile = "~{default_location}/~{outputfile}"
        File htmlfile = "~{default_location}/~{outputhtml}"
        File textfile = "~{default_location}/~{outputtext}"
        File xhtml = "~{default_location}/~{outputhtml}x"
    }
}

task normalize {
    input {
        File wigfile
        File xlsfile
        Boolean control = false
        String default_location = "Coverage_files"

        String outputfile = sub(basename(wigfile),'\.wig\.gz', '.RPM.wig')

        Int memory_gb = 5
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<

        mkdir -p ~{default_location} && cd ~{default_location}

        gunzip -c ~{wigfile} > ~{basename(wigfile,'.gz')}
        ln -s ~{basename(wigfile,'.gz')} thewig.wig
        ln -s ~{xlsfile} xlsfile.xls
        
        python <<CODE
        import subprocess
        
        if ("~{control}" == "true") :
                command1 = "grep 'tags after filtering in control' xlsfile.xls"
                command2 = "grep 'total tags in control' xlsfile.xls"
        else :
            command1 = "grep 'tags after filtering in treatment' xlsfile.xls"
            command2 = "grep 'total tags in treatment' xlsfile.xls"
        
        try:
            mappedreads = int(str(subprocess.check_output(command1,shell=True).strip()).split(': ')[1].split("'")[0].strip())
        except:
            mappedreads = 0
        if mappedreads <= 0:   
            mappedreads = int(str(subprocess.check_output(command2,shell=True).strip()).split(': ')[1].split("'")[0].strip()) 
        
        mappedreads = mappedreads/1000000

        print(mappedreads)
        
        inputwig = open("thewig.wig", 'r')
        Lines = inputwig.readlines()
        file1 = open("output.out", 'w') 
        
        for line in Lines :
            if line.startswith('track') or line.startswith('variable'):
                file1.write("%s" %(line))
            else:
                lines = line.split("\t")
                height = int(lines[1])/float(mappedreads)
                file1.write("%s\t%s\n" %(lines[0],height))
                
        CODE
        mv output.out ~{outputfile}
        gzip ~{outputfile}

    >>> 
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'abralab/seaseq:v2.0.0'
        cpu: ncpu
    }
    output {
        File norm_wig = "~{default_location}/~{outputfile}.gz"
    }
}


task peaksanno {
    input {
        File bedfile
	    File? summitfile
        String default_location = "PEAKSAnnotation"

        File gtffile
        File chromsizes

        Int memory_gb = 5
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<

        mkdir -p ~{default_location}

        cd ~{default_location}

        if [[ "~{gtffile}" == *"gz" ]]; then
            gunzip -c ~{gtffile} > ~{sub(basename(gtffile),'.gz','')}
        else
           ln -s ~{gtffile} ~{sub(basename(gtffile),'.gz','')}
        fi

        peaksanno.py \
        -p ~{bedfile} \
        ~{if defined(summitfile) then "-s " + summitfile else ""} \
        -g ~{sub(basename(gtffile),'.gz','')} \
        -c ~{chromsizes}

    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'abralab/seaseq:v2.0.0'
        cpu: ncpu
    }
    output {
        File? peak_promoters = "~{default_location}/peaks_within_promoter.regions.txt"
        File? peak_genebody = "~{default_location}/peaks_within_genebody.regions.txt"
        File? peak_window = "~{default_location}/peaks_within_window.regions.txt"
        File? peak_closest = "~{default_location}/centerofpeaks_closest.regions.txt"
        File? peak_comparison = "~{default_location}/peaks_compared_regions.peaks.txt"
        File? gene_comparison = "~{default_location}/peaks_compared_regions.genes.txt"
        File? pdf_comparison = "~{default_location}/peaks_compared_regions.distribution.pdf"

    }
}

task mergehtml {
    input {
        Array[File] htmlfiles
        Array[File] txtfiles
        String fastq_type = "Sample FASTQs"
        String default_location = "QC_files"
        String outputfile 

        Int memory_gb = 10
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<
        mkdir -p ~{default_location} && cd ~{default_location}

        #extract header information
        tail -n 101 /usr/local/bin/scripts/seaseq_overall.header > ~{outputfile}

        mergeoutput=$(cat ~{sep='; tail -n 1 ' htmlfiles})
        echo $mergeoutput >> ~{outputfile}
        sed -i "s/SEAseq ~{fastq_type} Report/SEAseq ~{fastq_type} Report/" ~{outputfile}

        echo $mergeoutput > ~{outputfile}x

        head -n 1 ~{txtfiles[0]} > ~{outputfile}.txt
        mergeoutput=$(tail -n 1 ~{sep='; echo "xxx"; tail -n 1 ' txtfiles})
        echo $mergeoutput | awk -F" xxx " '{for (i=1;i<=NF;i++) print $i}' >> ~{outputfile}.txt
        perl -pi -e 's/ /\t/g' ~{outputfile}.txt

    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'abralab/seaseq:v2.0.0'
        cpu: ncpu
    }
    output {
        File mergefile = "~{default_location}/~{outputfile}"
        File mergetxt = "~{default_location}/~{outputfile}.txt"
        File xhtml = "~{default_location}/~{outputfile}x"
    }
}

task linkname {
    input {
        String prefix
        File in_fastq

        Int memory_gb = 10
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<
        mv ~{in_fastq} ~{prefix}.fastq.gz
        
    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'abralab/seaseq:v2.0.0'
        cpu: ncpu
    }
    output {
        File out_fastq = "~{prefix}.fastq.gz"
    }
}
