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
        File? sampleqc_txt
        File? sampleqc_html
        File overallqc_txt
        File overallqc_html

        String outputfile = sub(basename(overallqc_html), 'stats.htmlx', 'seaseq_report.html')
        String outputtxt = sub(basename(overallqc_html), 'stats.htmlx', 'seaseq_report.txt')
        
        Int memory_gb = 5
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<

        # Printing header
        head -n 121 /usr/local/bin/scripts/seaseq_overall.header > ~{outputfile}
        if [ -f "~{sampleqc_html}" ]; then
            # Printing Sample Quality Reports
            echo '<h2>Sample FASTQs Quality Results</h2>' >> ~{outputfile}
            cat ~{sampleqc_html} >> ~{outputfile}
            echo -e '</table></div>\n<div class="body">' >> ~{outputfile}

            echo -e 'SEAseq Report\nSEAseq Quality Statistics and Evaluation Report\n\nSample FASTQs Quality Results' > ~{outputtxt}
            cat ~{sampleqc_txt} >> ~{outputtxt}
            echo -e '\n' >> ~{outputtxt}
        fi

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
        echo '</table>' >> ~{outputfile}
        echo -e 'Overall Quality Evaluation and Statistics Results' >> ~{outputtxt}
        cat ~{overallqc_txt} >> ~{outputtxt}
        if [ -f "~{sampleqc_html}" ]; then
            echo "<p><b>*</b> Peaks identified after Input/Control correction.</p>" >> ~{outputfile}
        fi
        echo '</div>' >> ~{outputfile}
        tail -n 13 /usr/local/bin/scripts/seaseq_overall.header >> ~{outputfile}
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
        File? bambed
        File? sppfile
        File fastqczip
        File? bamflag
        File? rmdupflag
        File? bkflag
        File? fastqmetrics
        File? countsfile
        File? peaksxls
        File? enhancers
        File? superenhancers

        String fastq_type = "Sample FASTQs"
        String default_location = "QC_files/STATS"
        String outputfile = sub(basename(fastqczip),'_fastqc.zip', '-stats.csv')
        String outputhtml = sub(basename(fastqczip),'_fastqc.zip', '-stats.html')
        String outputtext = sub(basename(fastqczip),'_fastqc.zip', '-stats.txt')
        String configml = sub(basename(fastqczip),'_fastqc.zip', '-config.ml')

        Int memory_gb = 10
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<

        mkdir -p ~{default_location}
        cd ~{default_location}

        evaluation-statistics.pl \
            -fqc ~{fastqczip} \
            ~{if defined(bambed) then "-b " + bambed else ""} \
            ~{if defined(sppfile) then "-s " + sppfile else ""} \
            ~{if defined(countsfile) then "-c " + countsfile else ""} \
            ~{if defined(peaksxls) then "-px " + peaksxls else ""} \
            ~{if defined(bamflag) then "-bamflag " + bamflag else ""} \
            ~{if defined(bkflag) then "-bkflag " + bkflag else ""} \
            ~{if defined(rmdupflag) then "-rmdupflag " + rmdupflag else ""} \
            ~{if defined(fastqmetrics) then "-fx " + fastqmetrics else ""} \
            ~{if defined(enhancers) then "-re " + enhancers else ""} \
            ~{if defined(superenhancers) then "-rs " + superenhancers else ""} \
            -outfile ~{outputfile}

        head -n 245 /usr/local/bin/scripts/seaseq_overall.header  | tail -n 123 > ~{outputhtml}
        cat ~{outputhtml}x >> ~{outputhtml}
        sed -i "s/SEAseq Sample FASTQ Report/SEAseq ~{fastq_type} Report/" ~{outputhtml}
        echo '</table></div>' >> ~{outputhtml}
        tail -n 13 /usr/local/bin/scripts/seaseq_overall.header >> ~{outputhtml}


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
        File configfile = "~{default_location}/~{configml}"
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
        continueOnReturnCode: [0, 1]
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
        head -n 245 /usr/local/bin/scripts/seaseq_overall.header  | tail -n 123 > ~{outputfile}

        mergeoutput=$(cat ~{sep='; tail -n 1 ' htmlfiles})
        echo $mergeoutput >> ~{outputfile}
        sed -i "s/SEAseq Sample FASTQ Report/SEAseq ~{fastq_type} Report/" ~{outputfile}
        echo '</table></div>' >> ~{outputfile}
        tail -n 13 /usr/local/bin/scripts/seaseq_overall.header >> ~{outputfile}

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


task concatstats {
    input {
        File sample_config
        File control_config
        File overall_config
        String outputfile = sub(basename(overall_config),'-config.ml','')
        String default_location = "QC_files"
        String fastq_type = "Sample FASTQs"

        Int memory_gb = 10
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<

        mkdir -p ~{default_location}
        cd ~{default_location}

        python <<CODE
        
        sample_content = open("~{sample_config}", 'r')
        control_content = open("~{control_config}", 'r')
        overall_content = open("~{overall_config}", 'r')
        htmlfile = "~{outputfile}-stats.htmlx"
        statsfile = htmlfile.replace("stats.htmlx", "stats.csv")
        textfile = htmlfile.replace("stats.htmlx", "stats.txt")

        writestatsfile = open(statsfile, 'w')

        #all statistics
        stats = {
            1 : 'Overall_Quality', 2 : 'Raw_Reads', 3 : 'Base_Quality', 
            4 : 'Sequence_Diversity', 5 : 'Aligned_Percent', 6 : 'Estimated_Fragment_Width', 
            7 : 'Estimated_Tag_Length', 8 : 'NRF', 9 : 'PBC', 10 : 'NSC', 11 : 'RSC', 
            12 : 'FRiP', 13 : 'Total_Peaks', 14 : 'Normalized_Peaks*', 
            15 : 'Linear_Stitched_Peaks', 16 : 'SE-like_Enriched_Regions'
        }

        #color scheme
        color = {
            -2 : "#FF0000", -1 : "#FF8C00", 0 : "#FFFF00", 1 : "#ADFF2F", 2 : "#008000"
        }

        #merged statistics
        mergestats = [1, 14, 15, 16]

        #initialize dictionaries
        SQCvalue = {}
        CQCvalue = {}
        OQCvalue = {}
        SQCscore = {}
        CQCscore = {}
        OQCscore = {}
        
        for line in sample_content:
            eachdata = line.rstrip('\n').split('\t')
            SQCvalue[eachdata[0]] = eachdata[1]
            SQCscore[eachdata[0]] = int(eachdata[2])
        for line in control_content:
            eachdata = line.rstrip('\n').split('\t')
            CQCvalue[eachdata[0]] = eachdata[1]
            CQCscore[eachdata[0]] = int(eachdata[2])
        for line in overall_content:
            eachdata = line.rstrip('\n').split('\t')
            OQCvalue[eachdata[0]] = eachdata[1]
            OQCscore[eachdata[0]] = int(eachdata[2])

        htmlheader = "<table class='results'><tr><th>DATA"
        textheader = "DATA"
        samplehtmlvalues = "<tr><td><center>SAMPLE</center></td>"
        controlhtmlvalues = "<tr><td><center>CONTROL</center></td>"
        sampletextvalues = "SAMPLE"
        controltextvalues = "CONTROL"

        for key in sorted(stats.keys()):
            #change space to underscore for txt fileprint(key)
            convertheader = stats[key].replace('_', ' ')
            textheader += "\t" + stats[key]
            htmlheader += "</th><th>" + convertheader
            if key in mergestats:
                if key == 14:
                    samplehtmlvalues += "<td rowspan='2' bgcolor='" + \
                                        color[OQCscore['Total_Peaks']] + "'><center>" + \
                                        OQCvalue['Total_Peaks'] + "</center></td>"
                    sampletextvalues += "\t" + OQCvalue['Total_Peaks']
                    controltextvalues += "\t"
                    writestatsfile.write(convertheader + "," + OQCvalue['Total_Peaks'] + "\n")
                else:
                    samplehtmlvalues += "<td rowspan='2' bgcolor='" + \
                                        color[OQCscore[stats[key]]] + "'><center>" + \
                                        OQCvalue[stats[key]] + "</center></td>"
                    sampletextvalues += "\t" + OQCvalue[stats[key]]
                    controltextvalues += "\t"
                    writestatsfile.write(convertheader + "," + OQCvalue[stats[key]] + "\n")
            else:
                samplehtmlvalues += "<td bgcolor='" + color[SQCscore[stats[key]]] + \
                                    "'><center>" + SQCvalue[stats[key]] + "</center></td>"
                controlhtmlvalues += "<td bgcolor='" + color[CQCscore[stats[key]]] + \
                                    "'><center>" + CQCvalue[stats[key]] + "</center></td>"
                sampletextvalues += "\t" + SQCvalue[stats[key]]
                controltextvalues += "\t" + CQCvalue[stats[key]]
                writestatsfile.write("Sample : " + convertheader + "," + SQCvalue[stats[key]] + "\n")
                writestatsfile.write("Control : " + convertheader + "," + SQCvalue[stats[key]] + "\n")

        htmlheader += "</th></tr>"
        samplehtmlvalues += "</tr>"
        controlhtmlvalues += "</tr>"

        writehtmlfile = open(htmlfile, 'w')
        writetextfile = open(textfile, 'w')

        writehtmlfile.write(htmlheader + "\n" + samplehtmlvalues + "\n" + controlhtmlvalues)
        writetextfile.write(textheader + "\n" + sampletextvalues + "\n" + controltextvalues + "\n")
        writestatsfile.write("\n" + '* Peaks identified after Input/Control correction')
        writetextfile.write("\n" + '* Peaks identified after Input/Control correction')

        CODE

        head -n 245 /usr/local/bin/scripts/seaseq_overall.header | tail -n 123 > ~{outputfile}-stats.html
        cat ~{outputfile}-stats.htmlx >> ~{outputfile}-stats.html
        echo "</table><p><b>*</b> Peaks identified after Input/Control correction.</p></div>" >> ~{outputfile}-stats.html
        tail -n 13 /usr/local/bin/scripts/seaseq_overall.header >> ~{outputfile}-stats.html
        sed -i "s/SEAseq Sample FASTQ Report/SEAseq Comprehensive Report/" ~{outputfile}-stats.html

    >>> 
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'abralab/seaseq:v2.0.0'
        cpu: ncpu
    }
    output {
        File statsfile = "~{default_location}/~{outputfile}-stats.csv"
        File htmlfile = "~{default_location}/~{outputfile}-stats.html"
        File textfile = "~{default_location}/~{outputfile}-stats.txt"
        File xhtml = "~{default_location}/~{outputfile}-stats.htmlx"
    }
}


task addreadme {
    input {
	    String default_location = "PEAKS_files"
        String output_file = "readme_peaks.txt"

        Int memory_gb = 2
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<
        mkdir -p ~{default_location}
        cd ~{default_location}
        echo 'SEAseq performs three peak calling algorithms, and they will be saved into the following descriptive folders:' > ~{output_file}
        echo -e '1.  NARROW_peaks   : For shorter or narrow regions of enrichment using MACS.' >> ~{output_file}
        echo '                     SEAseq performs three different peak calls:' >> ~{output_file}
        echo '                     a.  <samplename>-p9_kd-auto :  Peaks identified excluding duplicate tags.' >> ~{output_file}
        echo '                     b.  <samplename>-p9_kd-all  :  Peaks identified using duplicates to estimate signal.' >> ~{output_file}
        echo '                                                        This will be used to identify stitched peaks.' >> ~{output_file}
        echo '                     c.  <samplename>-nm         :  Peaks identified using a defined shift size (shiftsize=200).' >> ~{output_file}
        echo '                                                        This is used to generate an unbiased signal coverage plot.' >> ~{output_file}
        echo -e '\n2.  BROAD_peaks    : For broad peaks or broad domains using SICER.' >> ~{output_file}
        echo -e '\n3.  STITCHED_peaks : For clusters of stitched peaks identified using the ROSE program.\n' >> ~{output_file}

    >>>
    runtime {
	memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'abralab/seaseq:v2.0.0'
        cpu: ncpu
    }
    output {
	File readme_peaks = "~{default_location}/~{output_file}"
    }
}
