version 1.0
# PEAseq create fragments from BAM
# PEAseq BAM to BED

task fraggraph {
    input {
        File bamfile
        File chromsizes
        String bigwig = sub(basename(bamfile),"\.bam$", "\.w50\.FPM\.bw")
        String wig = sub(basename(bamfile),"\.bam$", "\.w50\.FPM\.wig")
        String tdf = sub(basename(bamfile),"\.bam$", "\.w50\.FPM\.tdf")
        String fragbam = sub(basename(bamfile),"\.bam$", "\.w50\.frag\.bam")
        String default_location = "BAM_files"

        Int memory_gb = 200
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<
        mkdir -p ~{default_location} && cd ~{default_location}
        ln -s ~{chromsizes} genome.chrom.sizes

        #namesort
        samtools sort -n \
            ~{bamfile} \
            -o ~{sub(basename(bamfile),"\.bam$", "\.ns\.bam")}

        #bamtobed
        bamToBed -bedpe \
            -i ~{sub(basename(bamfile),"\.bam$", "\.ns\.bam")} \
            > ~{sub(basename(bamfile),"\.bam$", "\.ns\.bed")}

        #sortbed
        awk -F\\t '{ 
            if($1==$4 && $1~"chr" && $6-$2<1000) 
            { print $1 "\t" $2 "\t" $6 "\t" $7 } 
            }' ~{sub(basename(bamfile),"\.bam$", "\.ns\.bed")} > ~{sub(basename(bamfile),"\.bam$", "\.ns\.sorted.bed")}

        #makewindows
        bedtools makewindows -w 50 \
            -g ~{chromsizes} \
            | sort -k1,1 -k2,2n \
            > ~{sub(basename(chromsizes),"\.tab", "-50bpwindows\.bed")}

        #bedtobam
        bedToBam \
            -i ~{sub(basename(bamfile),"\.bam$", "\.ns\.sorted.bed")} \
            -g ~{sub(basename(chromsizes),"\.tab", "-50bpwindows\.bed")} \
            > ~{sub(basename(bamfile),"\.bam$", "\.w50\.bam")}

        #position sort
        samtools sort \
            ~{sub(basename(bamfile),"\.bam$", "\.w50\.bam")} \
            -o ~{fragbam}

        #view fragments
        export mappedPE=$(samtools view -c -F 4 ~{fragbam})

        #create bedgraph
        intersectBed -c \
            -a ~{sub(basename(chromsizes),"\.tab", "-50bpwindows\.bed")} \
            -b ~{fragbam} \
            | awk -F'[\t]' -v mapped=$mappedPE '{print $1 "\t" $2 "\t" $3 "\t" $4*1000000/mapped} ' > ~{sub(basename(bamfile),"\.bam$", "\.w50\.FPM\.graph")}
 
        #create bigwig and tdf
        bedGraphToBigWig \
            ~{sub(basename(bamfile),"\.bam$", "\.w50\.FPM\.graph")} \
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
    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'abralab/peaseq_frag:latest'
        cpu: ncpu
    }
    output {
        File bigwigfile = "~{default_location}/~{bigwig}"
        File tdffile = "~{default_location}/~{tdf}"
        File wigfile = "~{default_location}/~{wig}"
        File fragbamfile = "~{default_location}/~{fragbam}"
    }
}

task pe_bamtobed {
    input {
        File bamfile
        String outputfile = sub(basename(bamfile),"\.bam$", "\.ns\.sorted.bed")

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
            -o ~{sub(basename(bamfile),"\.bam$", "\.ns\.bam")}

        #bamtobed
        bamToBed -bedpe \
            -i ~{sub(basename(bamfile),"\.bam$", "\.ns\.bam")} \
            > ~{sub(basename(bamfile),"\.bam$", "\.ns\.bed")}

        #sortbed
        awk -F\\t '{ 
            if($1==$4 && $1~"chr" && $6-$2<1000) 
            { print $1 "\t" $2 "\t" $6 "\t" $7 "\t255\t+" } 
            }' ~{sub(basename(bamfile),"\.bam$", "\.ns\.bed")} > ~{outputfile}

    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'abralab/peaseq_frag:latest'
        cpu: ncpu
    }
    output {
        File bedfile = "~{default_location}/~{outputfile}"
    }
}

