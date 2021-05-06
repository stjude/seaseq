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
        #================================================================================

        # 
        # GENERATING UCSC REFSEQ FILE
        #

        python <<CODE

        import re
        import sys
        
        gtf_input = open("~{sub(basename(gtffile),'.gz','')}",'r')
        gtf_name = "~{sub(basename(gtffile),'.gz','')}"
        print (gtf_input)
        refseq_output = open("genome_refseq.ucsc",'w')

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
            sys.exit("ERROR :\tGTF/GFF with either transcript/gene annotation is needed")
        print("NOTE :\tFeature type used is '%s'" %(feature))

        #organize gtf file
        gtf_input = open("~{sub(basename(gtffile),'.gz','')}",'r')
        for line in gtf_input:
            if not line.startswith('#'):
                lines = line.rstrip("\n").split("\t")
                newline = lines[8].split(';')
                if gtf_name.split('.')[-1] == 'gff' or gtf_name.split('.')[-1] == 'gff3':
                    if re.search(";transcript_id", lines[8]):
                        transcript_id = [s for s in newline if "transcript_id=" in s][0]
                    else:
                        transcript_id = newline[0]
                    if re.search(";gene_name", lines[8]):
                        gene_name = [s for s in newline if "gene_name=" in s][0]
                    else:
                        gene_name = newline[0]
                    transcript_id = re.sub('[\"\;]', '', transcript_id.split('=')[-1]) #clean transcript_id
                    gene_name = re.sub('[\"\;]', '', gene_name.split('=')[-1]) #clean gene_name
                elif gtf_name.split('.')[-1] == 'gtf':
                    if re.search("; transcript_id", lines[8]):
                        transcript_id = [s for s in newline if "transcript_id " in s][0]
                    else:
                        transcript_id = newline[0]
                    if re.search("; gene_name", lines[8]):
                        gene_name = [s for s in newline if "gene_name " in s][0]
                    else:
                        gene_name = newline[0]
                    transcript_id = re.sub('[\"\;]', '', transcript_id.split(' ')[-1]) #clean transcript_id
                    gene_name = re.sub('[\"\;]', '', gene_name.split(' ')[-1]) #clean gene_name
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
        ROSE_main.py -s $STITCH -t $TSS --custom genome_refseq.ucsc -i unionpeaks.gff -r $BAMFILE -o $OUTPUTDIR

        mkdir -p ~{default_location}
        mv $OUTPUTDIR/* ~{default_location}
    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'abralab/rose:v1.3.0'
        cpu: ncpu
    }
    output {
        File pngfile = "~{default_location}/unionpeaks_Plot_points.png"
        File? mapped_union = "~{default_location}/mappedGFF/unionpeaks_~{outputname}_MAPPED.gff"
        File? mapped_stitch = "~{default_location}/mappedGFF/unionpeaks_12.5KB_STITCHED_TSS_DISTAL_~{outputname}_MAPPED.gff"
        File enhancers = "~{default_location}/unionpeaks_AllStitched.table.txt"
        File super_enhancers = "~{default_location}/unionpeaks_SuperStitched.table.txt"
        File? gff_file = "~{default_location}/gff/unionpeaks.gff"
        File? gff_union = "~{default_location}/gff/unionpeaks_12.5KB_STITCHED_TSS_DISTAL.gff"
        File? union_enhancers = "~{default_location}/unionpeaks_Stitched_withSuper.bed"
        File? stitch_enhancers = "~{default_location}/unionpeaks_12.5KB_STITCHED_TSS_DISTAL_REGION_MAP.txt"
        File? e_to_g_enhancers = "~{default_location}/unionpeaks_AllStitched_REGION_TO_GENE.txt"
        File? g_to_e_enhancers = "~{default_location}/unionpeaks_AllStitched_GENE_TO_REGION.txt"
        File? e_to_g_super_enhancers = "~{default_location}/unionpeaks_SuperStitched_REGION_TO_GENE.txt"
        File? g_to_e_super_enhancers = "~{default_location}/unionpeaks_SuperStitched_GENE_TO_REGION.txt"
    }
}
