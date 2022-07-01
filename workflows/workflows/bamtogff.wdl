version 1.0

workflow bamtogff {
    input {
        File bamfile
        File bamindex
        File? control_bamfile
        File? control_bamindex

        File gtffile
        File chromsizes
        Int distance = 2000 #distance from site

        String samplename = basename(bamfile,'.bam')
        String default_location = "BAMDensity_files"
    }

    call bamtogff_gtftogenes {
        input :
            gtffile=gtffile,
            chromsizes=chromsizes,
            distance=distance
    }   

    call bamtogff_main as s_promoters {
        input :
            bamfile=bamfile,
            bamindex=bamindex,
            matrix_bins=100,
            annotation=bamtogff_gtftogenes.promoters,
            matrix_name="samplematrix-promoters.txt"
    }

    call bamtogff_main as s_genebody {
        input :
            bamfile=bamfile,
            bamindex=bamindex,
            matrix_bins=100,
            annotation=bamtogff_gtftogenes.genes,
            matrix_name="samplematrix-genebody.txt"
    }

    call bamtogff_main as s_upstream {
        input :
            bamfile=bamfile,
            bamindex=bamindex,
            matrix_bins=50,
            annotation=bamtogff_gtftogenes.upstream,
            matrix_name="samplematrix-upstream.txt"
    }

    call bamtogff_main as s_downstream {
        input :
            bamfile=bamfile,
            bamindex=bamindex,
            matrix_bins=50,
            annotation=bamtogff_gtftogenes.downstream,
            matrix_name="samplematrix-downstream.txt"
    }

    if (defined(control_bamfile)) {
        String string_controlbam = ""
        File control_bamfile_m = select_first([control_bamfile,string_controlbam])
        File control_bamindex_m = select_first([control_bamindex,string_controlbam])

        call bamtogff_main as c_promoters {
            input :
                bamfile=control_bamfile_m,
                bamindex=control_bamindex_m,
                matrix_bins=100,
                annotation=bamtogff_gtftogenes.promoters,
                matrix_name="inputmatrix-promoters.txt"
        }

        call bamtogff_main as c_genebody {
            input :
                bamfile=control_bamfile_m,
                bamindex=control_bamindex_m,
                matrix_bins=100,
                annotation=bamtogff_gtftogenes.genes,
                matrix_name="inputmatrix-genebody.txt"
        }

        call bamtogff_main as c_upstream {
            input :
                bamfile=control_bamfile_m,
                bamindex=control_bamindex_m,
                matrix_bins=50,
                annotation=bamtogff_gtftogenes.upstream,
                matrix_name="inputmatrix-upstream.txt"
        }

        call bamtogff_main as c_downstream {
            input :
                bamfile=control_bamfile_m,
                bamindex=control_bamindex_m,
                matrix_bins=50,
                annotation=bamtogff_gtftogenes.downstream,
                matrix_name="inputmatrix-downstream.txt"
        }
    }

    call bamtogff_plot {
        input :
            bamfile=bamfile,
            control_bamfile=control_bamfile,
            distance=distance,
            s_promoters=s_promoters.matrix_file,
            s_genebody=s_genebody.matrix_file,
            s_upstream=s_upstream.matrix_file,
            s_downstream=s_downstream.matrix_file,
            c_promoters=c_promoters.matrix_file,
            c_genebody=c_genebody.matrix_file,
            c_upstream=c_upstream.matrix_file,
            c_downstream=c_downstream.matrix_file,
            samplename=samplename,
            default_location=default_location
    }

    output {
        File s_matrices = bamtogff_plot.s_matrices
        File? c_matrices = bamtogff_plot.c_matrices
        File? densityplot = bamtogff_plot.densityplot
        File? pdf_gene = bamtogff_plot.pdf_gene
        File? pdf_h_gene = bamtogff_plot.pdf_h_gene
        File? png_h_gene = bamtogff_plot.png_h_gene
        File? jpg_h_gene = bamtogff_plot.jpg_h_gene
        File? pdf_promoters = bamtogff_plot.pdf_promoters
        File? pdf_h_promoters = bamtogff_plot.pdf_h_promoters
        File? png_h_promoters = bamtogff_plot.png_h_promoters
        File? jpg_h_promoters = bamtogff_plot.jpg_h_promoters
    }
}

task bamtogff_gtftogenes {
    
    input {
        File gtffile
        File chromsizes
        Int distance

        Int memory_gb = 5
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<
        if [[ "~{gtffile}" == *"gz" ]]; then
            gunzip -c ~{gtffile} > ~{sub(basename(gtffile),'.gz','')}
        else
           ln -s ~{gtffile} ~{sub(basename(gtffile),'.gz','')}
        fi

        BAM2GFF_gtftogenes.py -g ~{sub(basename(gtffile),'.gz','')} -c ~{chromsizes} -d ~{distance}
    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'ghcr.io/stjude/abralab/bamtogff:v1.2.2'
        cpu: ncpu
    }
    output {
        File promoters = "annotation/promoters.gff"
        File genes = "annotation/genes.gff"
        File upstream = "annotation/upstream.gff"
        File downstream = "annotation/downstream.gff"
    }

}

task bamtogff_main {

    input {
        File bamfile
        File bamindex
        File annotation
        Int matrix_bins

        String matrix_name = sub(basename(annotation),'.gff', '.txt')

        Int memory_gb = 5
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<
        ln -s ~{bamfile} ~{basename(bamfile)}
        ln -s ~{bamindex} ~{basename(bamindex)}
        
        BAM2GFF_main.py -b ~{basename(bamfile)} -i ~{annotation} -m ~{matrix_bins} -o ~{matrix_name}
    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'ghcr.io/stjude/abralab/bamtogff:v1.2.2'
        cpu: ncpu
    }
    output {
        File matrix_file = "~{matrix_name}"
    }

}

task bamtogff_plot {

    input {
        File bamfile
        File? control_bamfile
        Int distance    
        File s_promoters
        File s_genebody
        File s_upstream
        File s_downstream
        File? c_promoters
        File? c_genebody
        File? c_upstream
        File? c_downstream
        String default_location = "BAMDensity_files"

        String samplename = basename(bamfile,'.bam')
        
        Int memory_gb = 10
        Int max_retries = 1
        Int ncpu = 1
    }
    command <<<

        #create a custom R script to recreate plots
        RSCRIPT="densityplots.R"
        head -n 443 /opt/BAM2GFF-1.2.2/bin/BAM2GFF_plots.R > $RSCRIPT
        echo 'if ("pdftools" %in% rownames(installed.packages()) == FALSE){ install.packages("pdftools") }' >> $RSCRIPT
        echo '' >> $RSCRIPT
        echo '#=========================================' >> $RSCRIPT
        echo '' >> $RSCRIPT
        echo 'folder = "'sample_matrixfiles.zip'"' >> $RSCRIPT
        ~{if defined(control_bamfile) then "echo \'opt <- list(z=TRUE, c=\"input_matrixfiles.zip\")\' >> $RSCRIPT" else "echo \'opt <- list(z=TRUE, c=NA)\' >> $RSCRIPT"}
        echo 'samplename = "'~{samplename}'"' >> $RSCRIPT
        echo 'unzipped_folder = "UNZIPPED"' >> $RSCRIPT
        echo 'distance = round(~{distance}/1000,1)' >> $RSCRIPT
        echo '' >> $RSCRIPT
        echo '#=========================================' >> $RSCRIPT
        echo '' >> $RSCRIPT
        tail -n 127 /opt/BAM2GFF-1.2.2/bin/BAM2GFF_plots.R | head -n -3 >> $RSCRIPT 

        #moving sample matrix files to sample_matrixfiles folder
        mkdir -p ~{default_location} sample_matrixfiles
        cp ~{s_promoters} ~{s_genebody} ~{s_upstream} ~{s_downstream} sample_matrixfiles/

        cd sample_matrixfiles; zip -9r ../sample_matrixfiles.zip *; cd ..

        #creating plots w/ or w/o input bam files if provided
        if [ -f "~{control_bamfile}" ]; then
            mkdir -p input_matrixfiles
            cp ~{c_promoters} ~{c_genebody} ~{c_upstream} ~{c_downstream} input_matrixfiles/

            echo "Making PLOTS for ~{basename(bamfile)} read densities against INPUT"
            BAM2GFF_plots.R -n ~{samplename} -d ~{distance} -f sample_matrixfiles -c input_matrixfiles

            cd input_matrixfiles; zip -9r ../input_matrixfiles.zip *; cd ..

            mv input_matrixfiles.zip ~{default_location}

        else
            echo "Making PLOTS for ~{basename(bamfile)} read densities"
            BAM2GFF_plots.R -n ~{samplename} -d ~{distance} -f sample_matrixfiles

        fi

        if [ ! -f "~{samplename}-heatmap.entiregene.png" ]; then
            echo "using PDFtoPPM"
            #convert pdf to png #not needed since using pdftools
            pdftoppm ~{samplename}-heatmap.entiregene.pdf ~{samplename}-heatmap.entiregene -png -singlefile -r 300
            pdftoppm ~{samplename}-heatmap.promoters.pdf ~{samplename}-heatmap.promoters -png -singlefile -r 300
        fi

        mv $RSCRIPT sample_matrixfiles.zip *png *pdf *jpg ~{default_location}
        
        echo "Done!"

    >>>
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'ghcr.io/stjude/abralab/bamtogff:v1.2.2'
        cpu: ncpu
    }
    output {
        File s_matrices = "~{default_location}/sample_matrixfiles.zip"
        File? c_matrices = "~{default_location}/input_matrixfiles.zip"
        File? densityplot = "~{default_location}/densityplots.R"
        File? pdf_gene = "~{default_location}/~{samplename}-entiregene.pdf"
        File? pdf_h_gene = "~{default_location}/~{samplename}-heatmap.entiregene.pdf"
        File? png_h_gene = "~{default_location}/~{samplename}-heatmap.entiregene.png"
        File? jpg_h_gene = "~{default_location}/~{samplename}-heatmap.entiregene.jpg"
        File? pdf_promoters = "~{default_location}/~{samplename}-promoters.pdf"
        File? pdf_h_promoters = "~{default_location}/~{samplename}-heatmap.promoters.pdf"
        File? png_h_promoters = "~{default_location}/~{samplename}-heatmap.promoters.png"
        File? jpg_h_promoters = "~{default_location}/~{samplename}-heatmap.promoters.jpg"
    }

}

