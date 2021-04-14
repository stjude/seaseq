version 1.0

task bamtogff {

    input {
        File bamfile
        File bamindex
        String default_location = "BAMDensity_files"

        File gtffile
        File chromsizes
        Int distance = 2000 #distance from site

        String samplename = basename(bamfile,'.bam')

        Int memory_gb = 5
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<
        mkdir -p ~{default_location}

        ln -s ~{bamfile} ~{basename(bamfile)}
        ln -s ~{bamindex} ~{basename(bamindex)}

        if [[ "~{gtffile}" == *"gz" ]]; then
            gunzip -c ~{gtffile} > ~{sub(basename(gtffile),'.gz','')}
        else
           ln -s ~{gtffile} ~{sub(basename(gtffile),'.gz','')}
        fi

        GTFFILE=~{sub(basename(gtffile),'.gz','')}
        BAMFILE=~{basename(bamfile)}
        CHROMSIZES=~{chromsizes}
        SAMPLENAME=~{samplename}
        MATRIXFILES="matrixfiles"
        RSCRIPT="densityplots.R"

        echo "#############################################"
        echo "######            BAM2GFF v1           ######"
        echo "#############################################"

        echo "BAM file: $BAMFILE"
        echo "Sample Name: $SAMPLENAME"

        echo "Extracting GENE regions"
        BAM2GFF_gtftogenes.py -g $GTFFILE -c $CHROMSIZES -d ~{distance}

        mkdir -p $MATRIXFILES
        echo "Working on Promoter Region"
        BAM2GFF_main.py -b $BAMFILE -i annotation/promoters.gff -m 100 -o $MATRIXFILES/promoters.txt

        echo "Working on GeneBody Region"
        BAM2GFF_main.py -b $BAMFILE -i annotation/genes.gff -m 100 -o $MATRIXFILES/genebody.txt

        echo "Working on Upstream Region"
        BAM2GFF_main.py -b $BAMFILE -i annotation/upstream.gff -m 50 -o $MATRIXFILES/upstream.txt

        echo "Working on Downstream Region"
        BAM2GFF_main.py -b $BAMFILE -i annotation/downstream.gff -m 50 -o $MATRIXFILES/downstream.txt

        echo "BAM2GFF_plots.R $SAMPLENAME"
        BAM2GFF_plots.R -n $SAMPLENAME -d ~{distance} -f $MATRIXFILES

        #create a custom R script to recreate plots
        head -n 461 /opt/BAM2GFF-1.2.0/bin/BAM2GFF_plots.R > $RSCRIPT
        echo 'folder = "'$MATRIXFILES'"' >> $RSCRIPT
        echo 'samplename = "'$SAMPLENAME'"' >> $RSCRIPT
        echo 'distance = round(~{distance}/1000,1)' >> $RSCRIPT
        tail -n 53 /opt/BAM2GFF-1.2.0/bin/BAM2GFF_plots.R >> $RSCRIPT 
        echo "Done!"

        mv $RSCRIPT $MATRIXFILES *png *pdf ~{default_location}
    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'abralab/bamtogff:v1.2.0'
        cpu: ncpu
    }
    output {
        File m_downstream = "~{default_location}/matrixfiles/downstream.txt"
	File m_upstream = "~{default_location}/matrixfiles/upstream.txt"
	File m_genebody = "~{default_location}/matrixfiles/genebody.txt"
	File m_promoters = "~{default_location}/matrixfiles/promoters.txt"
        File densityplot = "~{default_location}/densityplots.R"
        File? pdf_gene = "~{default_location}/~{samplename}-entiregene.pdf"
        File? pdf_h_gene = "~{default_location}/~{samplename}-heatmap.entiregene.pdf"
        File? png_h_gene = "~{default_location}/~{samplename}-heatmap.entiregene.png"
        File? pdf_promoters = "~{default_location}/~{samplename}-promoters.pdf"
        File? pdf_h_promoters = "~{default_location}/~{samplename}-heatmap.promoters.pdf"
        File? png_h_promoters = "~{default_location}/~{samplename}-heatmap.promoters.png"
    }
}
