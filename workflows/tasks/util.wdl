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
        docker: 'abralab/seaseq:v1.1.0'
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
        docker: 'abralab/seaseq:v1.1.0'
        cpu: ncpu
    }
    output {
        File flankbedfile = "~{default_location}/~{outputfile}"
    }
}

task summarystats {
    input {
        File? bambed
        File? sppfile
        File? countsfile
        File? peaksxls
        File? bamflag
        File? rmdupflag
        File? bkflag
        File fastqczip
        File? fastqmetrics
        File? enhancers
        File? superenhancers

        String default_location = "QC_files/STATS"
        String outputfile = sub(basename(fastqczip),'_fastqc.zip', '_seaseq-summary-stats.csv')
        String outputhtml = sub(basename(fastqczip),'_fastqc.zip', '_seaseq-summary-stats.html')
        String outputtext = sub(basename(fastqczip),'_fastqc.zip', '_seaseq-summary-stats.txt')

        Int memory_gb = 10
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<
        mkdir -p ~{default_location}

        cd ~{default_location}

        summaryfacts.pl \
            ~{if defined(bambed) then "-b " + bambed else ""} \
            ~{if defined(sppfile) then "-s " + sppfile else ""} \
            ~{if defined(countsfile) then "-c " + countsfile else ""} \
            ~{if defined(peaksxls) then "-px " + peaksxls else ""} \
            ~{if defined(bamflag) then "-bamflag " + bamflag else ""} \
            ~{if defined(bkflag) then "-bkflag " + bkflag else ""} \
            ~{if defined(fastqczip) then "-fqc " + fastqczip else ""} \
            ~{if defined(fastqmetrics) then "-fx " + fastqmetrics else ""} \
            ~{if defined(enhancers) then "-re " + enhancers else ""} \
            ~{if defined(superenhancers) then "-rs " + superenhancers else ""} \
            -outfile ~{outputfile}
    >>> 
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'abralab/seaseq:v1.1.0'
        cpu: ncpu
    }
    output {
        File statsfile = "~{default_location}/~{outputfile}"
        File htmlfile = "~{default_location}/~{outputhtml}"
        File textfile = "~{default_location}/~{outputtext}"
    }
}


task normalize {
    input {
        File wigfile
        File xlsfile
        String default_location = "PEAKDisplay_files"

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
        
        command = "grep 'tags after filtering in treatment' xlsfile.xls"
        try:
            mappedreads = int(str(subprocess.check_output(command,shell=True).strip()).split(': ')[1].split("'")[0].strip())
        except:
            mappedreads = 0
        if mappedreads <= 0:
            command = "grep 'total tags in treatment' xlsfile.xls"
            mappedreads = int(str(subprocess.check_output(command,shell=True).strip()).split(': ')[1].split("'")[0].strip())  

        mappedreads = mappedreads/1000000
        
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
        docker: 'abralab/seaseq:v1.1.0'
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
        docker: 'abralab/seaseq:v1.1.0'
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
        String default_location = "QC_files"
        String outputfile 

        Int memory_gb = 10
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<
        mkdir -p ~{default_location} && cd ~{default_location}

        mergeoutput=$(cat ~{sep='; tail -n 1 ' htmlfiles})
        echo $mergeoutput > ~{outputfile}
    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'abralab/seaseq:v1.1.0'
        cpu: ncpu
    }
    output {
        File mergefile = "~{default_location}/~{outputfile}"
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
        docker: 'abralab/seaseq:v1.1.0'
        cpu: ncpu
    }
    output {
        File out_fastq = "~{prefix}.fastq.gz"
    }
}
