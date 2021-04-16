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

        if [[ "~{gtffile}" == *"gz" ]]; then
            gunzip -c ~{gtffile} > ~{sub(basename(gtffile),'.gz','')}
        else
           ln -s ~{gtffile} ~{sub(basename(gtffile),'.gz','')}
        fi


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
        #================================================================================

        # 
        # GENERATING UCSC REFSEQ FILE
        #
        mkdir -p annotation

        python <<CODE

        import re
        import sys
        
        gtf_input = open("~{sub(basename(gtffile),'.gz','')}",'r')
        print (gtf_input)
        refseq_output = open("annotation/~{species}_refseq.ucsc",'w')

        refseq_output.write("#bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\tX\tX\tX\t\tX\tname2\n")

        #reading the feature type to be used
        feature_dict = {}
        for line in gtf_input:
            if not line.startswith('#'):
                lines = line.split("\t")
                feature_dict[lines[2]] = lines[2]

        if 'transcript' in feature_dict:
            feature = "transcript"
        elif 'gene' in feature_dict:
            feature = "gene"
        else:
            sys.exit("ERROR :\tGTF with either transcript/gene annotation is needed")
        print("NOTE :\tFeature type used is '%s'" %(feature))

        #organize gtf file
        gtf_input = open("~{sub(basename(gtffile),'.gz','')}",'r')
        for line in gtf_input:
            if not line.startswith('#'):
                lines = line.split("\t")
                newline = lines[8].split(' ')
                gene_name = re.sub('[\"\;]', '', newline[1]) #clean gene_name
                transcript_id = re.sub('[\"\;]', '', newline[3]) #clean transcript_id
                if lines[2] == feature:
                    if re.match("chr", lines[0]):
                        new_chr_name = lines[0]
                    else:
                        new_chr_name = "chr"+lines[0]
                    
                    results = ("0\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t.\t.\t.\t.\t{7}".format(transcript_id, new_chr_name, lines[6], lines[3], lines[4], lines[3], lines[4], gene_name))

                    #extract chromosomal coordinates and store in outputfiles
                    refseq_output.write(results + "\n")
        CODE

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
        docker: 'abralab/rose:v1.2.0'
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
