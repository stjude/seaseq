version 1.0

task rose {

    input {
        File bamfile
        File bamindex
        File bedfile_auto
        File bedfile_all
        File gtffile

        String outputdir = "ROSE_out"
        String default_location = "PEAKS_files/STITCHED_REGIONS"
        String outputname = basename(bamfile)

        String feature = "gene"
        String species = "hg19"
        Int tss = 2000
        Int stitch = 12500

        Int memory_gb = 20
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<
        ln -s ~{bamfile} ~{basename(bamfile)}
        ln -s ~{bamindex} ~{basename(bamindex)}

        GTFFILE=~{gtffile}
        FEATURE=~{feature}
        SPECIES=~{species}
        BAMFILE=~{basename(bamfile)}
        TSS=~{tss}
        STITCH=~{stitch}
        FILEA=~{bedfile_auto}
        FILEB=~{bedfile_all}
        OUTPUTDIR=~{outputdir}

        echo "#############################################"
        echo "######             ROSE v1             ######"
        echo "#############################################"

        echo "Input Bed File A: $FILEA"
        echo "Input Bed File B: $FILEB"
        echo "BAM file: $BAMFILE"
        echo "Output directory: $OUTPUTDIR"
        echo "Species: $SPECIES"
        echo "Feature type: $FEATURE"
        #================================================================================
        # 
        # GENERATING UCSC REFSEQ FILE
        #
        mkdir -p annotation
        echo -e "#bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\tX\tX\tX\t\tX\tname2" > annotation/$SPECIES"_refseq.ucsc"

        if [[ $FEATURE == "gene" ]]; then
        awk -F'[\t ]' '{
            if($3=="gene")
                print "0\t" $14 "\tchr" $1 "\t" $7 "\t" $4 "\t" $5 "\t" $4 "\t" $5 "\t.\t.\t.\t.\t" $18}' $GTFFILE | sed s/\"//g >> annotation/$SPECIES"_refseq.ucsc"

        elif [[ $FEATURE == "transcript" ]]; then
        awk -F'[\t ]' '{
            if($3=="transcript")
                print "0\t" $14 "\tchr" $1 "\t" $7 "\t" $4 "\t" $5 "\t" $4 "\t" $5 "\t.\t.\t.\t.\t" $18}' $GTFFILE | sed s/\"//g >> annotation/$SPECIES"_refseq.ucsc"
        fi
        echo "Annotation file: "$SPECIES"_refseq.ucsc"

        #
        # GENERATING MERGE BED FILES
        #
        cat $FILEA $FILEB | sort -k1,1 -k2,2n | mergeBed -i - | awk -F\\t '{print $1 "\t" NR "\t\t" $2 "\t" $3 "\t\t.\t\t" NR}' > unionpeaks.gff
        echo "Merge Bed file: unionpeaks.gff"
        echo

        #
        # ROSE CALLER
        #
        ROSE_main.py -s $STITCH -t $TSS -g $SPECIES -i unionpeaks.gff -r $BAMFILE -o $OUTPUTDIR

        mkdir -p ~{default_location}
        mv $OUTPUTDIR/* ~{default_location}
    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'madetunj/rose:v1.2.0'
        cpu: ncpu
    }
    output {
        File pngfile = "~{default_location}/unionpeaks_Plot_points.png"
        File? mapped_union = "~{default_location}/mappedGFF/unionpeaks_~{outputname}_MAPPED.gff"
        File? mapped_stitch = "~{default_location}/mappedGFF/unionpeaks_12.5KB_STITCHED_TSS_DISTAL_~{outputname}_MAPPED.gff"
        File enhancers = "~{default_location}/unionpeaks_AllEnhancers.table.txt"
        File super_enhancers = "~{default_location}/unionpeaks_SuperEnhancers.table.txt"
        File? gff_file = "~{default_location}/gff/unionpeaks.gff"
        File? gff_union = "~{default_location}/gff/unionpeaks_12.5KB_STITCHED_TSS_DISTAL.gff"
        File? union_enhancers = "~{default_location}/unionpeaks_Enhancers_withSuper.bed"
        File? stitch_enhancers = "~{default_location}/unionpeaks_12.5KB_STITCHED_TSS_DISTAL_ENHANCER_REGION_MAP.txt"
        File? e_to_g_enhancers = "~{default_location}/unionpeaks_AllEnhancers_ENHANCER_TO_GENE.txt"
        File? g_to_e_enhancers = "~{default_location}unionpeaks_AllEnhancers_GENE_TO_ENHANCER.txt"
        File? e_to_g_super_enhancers = "~{default_location}/unionpeaks_SuperEnhancers_ENHANCER_TO_GENE.txt"
        File? g_to_e_super_enhancers = "~{default_location}/unionpeaks_SuperEnhancers_GENE_TO_ENHANCER.txt"
    }
}
