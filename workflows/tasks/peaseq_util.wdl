version 1.0
# PEAseq create fragments from BAM
# PEAseq BAM to BED

task fraggraph {
    input {
        File bamfile
        File chromsizes
        String bigwig = sub(basename(bamfile),".bam$", ".w50.FPM.bw")
        String wig = sub(basename(bamfile),".bam$", ".w50.FPM.wig")
        String tdf = sub(basename(bamfile),".bam$", ".w50.FPM.tdf")
        String fragbam = sub(basename(bamfile),".bam$", ".frag.bam")
        String bampebed = sub(basename(bamfile),".bam$", ".bam2bedpe.bed")
        String fragsizes = sub(basename(bamfile),".bam$", ".fragments.png")
        String default_location = "BAM_files"
        String bam_location = "BAM_files"
        String annotation_location = "Annotation"

        Int memory_gb = 50
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<
        cwd=$(pwd)
        mkdir -p ~{default_location} ~{bam_location} ~{annotation_location}
        cd ~{default_location}
        ln -s ~{chromsizes} genome.chrom.sizes

        #namesort
        samtools sort -n \
            ~{bamfile} \
            -o ~{sub(basename(bamfile),".bam$", ".ns.bam")}

        #bamtobed
        bamToBed -bedpe \
            -i ~{sub(basename(bamfile),".bam$", ".ns.bam")} \
            > ~{sub(basename(bamfile),".bam$", ".ns2bedpe.bedpe")}

        #filterbedpe
        awk -F\\t '{
            if($1==$4 && $1~"chr")
            { print $1 "\t" $2 "\t" $6 "\t" $7 "\t255\t+" }
            }' ~{sub(basename(bamfile),".bam$", ".ns2bedpe.bedpe")} > ~{bampebed}

        #fragmentsgraph
        python3 <<CODE

        import sys
        import matplotlib.pyplot as plt
        import pandas as pd
        import seaborn as sns

        bed_input = open("~{bampebed}",'r')

        ALLSIZES = []
        for line in bed_input:
            line = line.strip('\n').split('\t')
            size = int(line[2]) - int(line[1])
            ALLSIZES.append(size)

        df = pd.DataFrame({'Value':ALLSIZES})

        medians = df["Value"].median()
        maximum = df["Value"].max()
        means = df["Value"].mean()
        quantile3 = float(df["Value"].quantile([.75]))
        minimum = df["Value"].min()
        quantile1 = float(df["Value"].quantile([.5]))
        nobs = "Q2 = " + str(round(medians,2))
        nobs_max = "maximum = " + str(maximum)
        nobs_mean = "mean = " + str(round(means,2))
        nobs_min = "minimum = " + str(minimum)
        nobs_q1 = "Q1 = " + str(round(quantile1,2))
        nobs_q3 = "Q3 = " + str(round(quantile3,2))

        bbox = dict(boxstyle ="round", fc ="0.8")
        arrowprops = dict(
            arrowstyle = "-", color = 'black',
            connectionstyle = "angle, angleA = 0, \
            angleB = 90, rad = 10")

        ax = sns.violinplot(x=df["Value"], palette="pastel")
        ax.annotate(nobs_q1, xy = (quantile1,0), xytext=(quantile1,-0.3), bbox=bbox, arrowprops=arrowprops, size='small', weight='semibold')
        ax.annotate(nobs, xy = (medians,0), xytext=(medians,-0.2), bbox=bbox, arrowprops=arrowprops, size='small', weight='semibold')
        ax.annotate(nobs_q3, xy = (quantile3,0), xytext=(quantile3,-0.1), bbox=bbox, arrowprops=arrowprops, size='small', weight='semibold')
        ax.annotate(nobs_min, xy = (minimum,0), xytext=(minimum,0.3), bbox=bbox, arrowprops=arrowprops, size='small', weight='semibold')
        ax.annotate(nobs_mean, xy = (means,0), xytext=(means,0.2), bbox=bbox, arrowprops=arrowprops, size='small', weight='semibold')
        ax.annotate(nobs_max, xy = (maximum,0), xytext=(maximum,0.1), bbox=bbox, arrowprops=arrowprops, size='small', weight='semibold', ha='right')

        ax.set_xlabel("Fragments sizes (bp)")
        ax.set_title("~{sub(basename(bamfile),".sorted.*$", "")}",size='larger',weight='bold')

        plt.tick_params(left = False)
        plt.savefig("~{fragsizes}")

        CODE

        #makewindows
        bedtools makewindows -w 50 \
            -g ~{chromsizes} \
            | sort -k1,1 -k2,2n \
            > ~{sub(basename(chromsizes),".tab", "-50bpwindows.bed")}

        #bedtobam
        bedToBam \
            -i ~{bampebed} \
            -g ~{chromsizes} \
            > ~{sub(basename(bamfile),".bam$", ".bam2bedpe.bam")}

        #position sort
        samtools sort \
            ~{sub(basename(bamfile),".bam$", ".bam2bedpe.bam")} \
            -o ~{fragbam}

        #view fragments
        export mappedPE=$(wc -l ~{bampebed} | awk -F" " '{print $1}')

        #create bedgraph
        intersectBed -c \
            -a ~{sub(basename(chromsizes),".tab", "-50bpwindows.bed")} \
            -b ~{bampebed} \
            | awk -F'[\t]' -v mapped=$mappedPE '{print $1 "\t" $2 "\t" $3 "\t" $4*1000000/mapped}' \
            > ~{sub(basename(bamfile),".bam$", ".w50.FPM.graph")}

        #create bigwig and tdf
        bedGraphToBigWig \
            ~{sub(basename(bamfile),".bam$", ".w50.FPM.graph")} \
            ~{chromsizes} \
            ~{bigwig}

        bigWigToWig \
            ~{bigwig} \
            ~{wig}

        igvtools \
            toTDF \
            ~{wig} \
            ~{tdf} \
            genome.chrom.sizes

        gzip ~{wig}

        # move fragbam to new location
        mv ~{fragbam} $cwd/~{bam_location}
        mv ~{fragsizes} $cwd/~{annotation_location}
    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'ghcr.io/stjude/seaseq/data_processing:v1.0.0'
        cpu: ncpu
    }
    output {
        File bigwigfile = "~{default_location}/~{bigwig}"
        File tdffile = "~{default_location}/~{tdf}"
        File wigfile = "~{default_location}/~{wig}.gz"
        File fragbamfile = "~{bam_location}/~{fragbam}"
        File bedpefile = "~{default_location}/~{bampebed}"
        File fragsizepng = "~{annotation_location}/~{fragsizes}"
    }
}

task pe_bamtobed {
    input {
        File bamfile
        String outputfile = basename(bamfile) + "2bedpe.bed"
        Boolean sicer = false

        String default_location = "BAM_files"

        Int memory_gb = 20
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<
        mkdir -p ~{default_location} && cd ~{default_location}

        #namesort
        samtools sort -n \
            ~{bamfile} \
            -o ~{sub(basename(bamfile),".bam$", ".ns.bam")}

        #bamtobed
        bamToBed -bedpe \
            -i ~{sub(basename(bamfile),".bam$", ".ns.bam")} \
            > ~{sub(basename(bamfile),".bam$", ".ns.bed")}

        #sortbed
        awk -F\\t '{
            if($1==$4 && $1~"chr")
            { print $1 "\t" $2 "\t" $6 "\t" $7 "\t255\t+" }
            }' ~{sub(basename(bamfile),".bam$", ".ns.bed")} > ~{outputfile}
        #fi

    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'ghcr.io/stjude/seaseq/data_processing:v1.0.0'
        cpu: ncpu
    }
    output {
        File bedfile = "~{default_location}/~{outputfile}"
    }
}

task pairedend_summaryreport {
    input {
        File? controlqc_se_txt
        File? controlqc_se_html
        File? sampleqc_se_txt
        File? sampleqc_se_html
        File overallqc_se_txt
        File overallqc_se_html
        File? controlqc_pe_txt
        File? controlqc_pe_html
        File? sampleqc_pe_txt
        File? sampleqc_pe_html
        File overallqc_pe_txt
        File overallqc_pe_html

        String outputfile
        String outputtxt = sub(outputfile, '.html', '.txt')

        Int memory_gb = 5
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<

        # Printing header
        head -n 121 /usr/local/bin/seaseq_overall.header > ~{outputfile}
        sed -i "s/SEAseq Report/PEAseq Report/" ~{outputfile}
        sed -i "s/SEAseq Quality/PEAseq Quality/" ~{outputfile}
        if [ -f "~{sampleqc_se_html}" ]; then
            # Printing Sample Quality Reports
            echo '<h2>Sample FASTQs Quality Results</h2>' >> ~{outputfile}
            echo '<h3>Individual FASTQs</h3>' >> ~{outputfile}
            cat ~{sampleqc_se_html} >> ~{outputfile}
            echo -e '</table>' >> ~{outputfile}

            echo -e 'PEAseq Report\nSEAseq Quality Statistics and Evaluation Report\n\nSample FASTQs Quality Results\nIndividual FASTQs' > ~{outputtxt}
            cat ~{sampleqc_se_txt} >> ~{outputtxt}

            if [ -f "~{sampleqc_pe_html}" ]; then
                echo '<br><h3>After Paired End Reference Mapping</h3>' >> ~{outputfile}
                cat ~{sampleqc_pe_html} >> ~{outputfile}
                echo -e '</table></div>\n<div class="body">' >> ~{outputfile}
                echo -e '\nAfter Paired End Reference Mapping' > ~{outputtxt}
                cat ~{sampleqc_pe_txt} >> ~{outputtxt}
                echo -e '\n' >> ~{outputtxt}
            else
                echo -e '</div>\n<div class="body">' >> ~{outputfile}
                echo -e '\n' >> ~{outputtxt}
            fi

        fi

        if [ -f "~{controlqc_se_html}" ]; then
            # Printing Control Quality Reports
            echo '<h2>Control FASTQs Quality Results</h2>' >> ~{outputfile}
            echo '<h3>Individual FASTQs</h3>' >> ~{outputfile}
            cat ~{controlqc_se_html} >> ~{outputfile}
            echo -e '</table>' >> ~{outputfile}

            echo 'Control FASTQs Quality Results\nIndividual FASTQs' >> ~{outputtxt}
            cat ~{controlqc_se_txt} >> ~{outputtxt}

            if [ -f "~{controlqc_pe_html}" ]; then
                echo '<br><h3>After Paired End Reference Mapping</h3>' >> ~{outputfile}
                cat ~{controlqc_pe_html} >> ~{outputfile}
                echo -e '</table></div>\n<div class="body">' >> ~{outputfile}
                echo -e '\nAfter Paired End Reference Mapping' > ~{outputtxt}
                cat ~{controlqc_pe_txt} >> ~{outputtxt}
                echo -e '\n' >> ~{outputtxt}
            else
                echo -e '</div>\n<div class="body">' >> ~{outputfile}
                echo -e '\n' >> ~{outputtxt}
            fi

        fi

        # Printing Overall Quality Reports
        echo '<h2>Overall Quality Evaluation and Statistics Results</h2><p>' >> ~{outputfile}
        echo '<h3>Single End Mode</h3>' >> ~{outputfile}
        cat ~{overallqc_se_html} >> ~{outputfile}
        echo '</table>' >> ~{outputfile}
        echo '<br><h3>Paired End Mode</h3>' >> ~{outputfile}
        cat ~{overallqc_pe_html} >> ~{outputfile}
        echo '</table>' >> ~{outputfile}

        echo -e 'Overall Quality Evaluation and Statistics Results\nSingle End mode' >> ~{outputtxt}
        cat ~{overallqc_se_txt} >> ~{outputtxt}
        echo -e '\nPaired End mode' >> ~{outputtxt}
        cat ~{overallqc_pe_txt} >> ~{outputtxt}

        if [ -f "~{controlqc_se_html}" ]; then
            echo "<p><b>*</b> Peaks identified after Input/Control correction.</p>" >> ~{outputfile}
        fi
        echo '</div>' >> ~{outputfile}
        tail -n 13 /usr/local/bin/seaseq_overall.header >> ~{outputfile}
        echo -e '\n' >> ~{outputtxt}

    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'ghcr.io/stjude/seaseq/scripts:v2.4.0'
        cpu: ncpu
    }
    output {
        File summaryhtml = "~{outputfile}"
        File summarytxt = "~{outputtxt}"
    }
}

task pe_mergehtml {
    input {
        Array[File] se_htmlfiles
        Array[File] se_txtfiles
        Array[File] pe_htmlfiles
        Array[File] pe_txtfiles
        String fastq_type = "PEAseq Sample FASTQs"
        String default_location = "QC_files"
        String outputfile 
        String outputtxt = sub(outputfile, '.html', '.txt')
        Boolean peaseq = false

        Int memory_gb = 10
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<
        mkdir -p ~{default_location} && cd ~{default_location}

        #extract header information
        head -n 245 /usr/local/bin/seaseq_overall.header  | tail -n 123 > ~{outputfile}
        sed -i "s/SEAseq Report/PEAseq Report/" ~{outputfile}
        sed -i "s/SEAseq Quality/PEAseq Quality/" ~{outputfile}

        echo '<h3>Individual FASTQs</h3>' >> ~{outputfile}
            
        se_mergeoutput=$(cat ~{sep='; tail -n 1 ' se_htmlfiles})
        echo $se_mergeoutput >> ~{outputfile}
        sed -i "s/SEAseq Sample FASTQ Report/~{fastq_type} Report/" ~{outputfile}
        echo '</table>' >> ~{outputfile}
        
        echo '<br><h3>After Paired End Reference Mapping</h3>' >> ~{outputfile}
        pe_mergeoutput=$(cat ~{sep='; tail -n 1 ' pe_htmlfiles})
        echo $pe_mergeoutput >> ~{outputfile}
        echo '</table></div>' >> ~{outputfile}
        
        tail -n 13 /usr/local/bin/seaseq_overall.header >> ~{outputfile}

        #working on TXTfile
        head -n 1 ~{se_txtfiles[0]} > ~{outputtxt}
        se_mergeoutput=$(tail -n 1 ~{sep='; echo "xxx"; tail -n 1 ' se_txtfiles})
        pe_mergeoutput=$(tail -n 1 ~{sep='; echo "xxx"; tail -n 1 ' pe_txtfiles})
        echo $se_mergeoutput | awk -F" xxx " '{for (i=1;i<=NF;i++) print $i}' >> ~{outputtxt}
        echo -e '\nAfter Paired End Reference Mapping\n' >> ~{outputtxt}
        echo $pe_mergeoutput | awk -F" xxx " '{for (i=1;i<=NF;i++) print $i}' >> ~{outputtxt}
        perl -pi -e 's/ /\t/g' ~{outputtxt}

    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'ghcr.io/stjude/seaseq/scripts:v2.4.0'
        cpu: ncpu
    }
    output {
        File mergefile = "~{default_location}/~{outputfile}"
        File mergetxt = "~{default_location}/~{outputtxt}"
    }
}
