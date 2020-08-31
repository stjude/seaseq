version 1.0
import "https://raw.githubusercontent.com/madetunj/wdl-seaseq/master/wdl/fastqc.wdl"
import "https://raw.githubusercontent.com/madetunj/wdl-seaseq/master/wdl/bedtools.wdl"
import "https://raw.githubusercontent.com/madetunj/wdl-seaseq/master/wdl/bowtie.wdl"
import "https://raw.githubusercontent.com/madetunj/wdl-seaseq/master/wdl/samtools.wdl"
import "https://raw.githubusercontent.com/madetunj/wdl-seaseq/master/wdl/macs.wdl" 
import "https://raw.githubusercontent.com/madetunj/wdl-seaseq/master/wdl/bamtogff.wdl"
import "https://raw.githubusercontent.com/madetunj/wdl-seaseq/master/wdl/sicer.wdl"
import "https://raw.githubusercontent.com/madetunj/wdl-seaseq/master/wdl/motifs.wdl"
import "https://raw.githubusercontent.com/madetunj/wdl-seaseq/master/wdl/rose.wdl"
import "https://raw.githubusercontent.com/madetunj/wdl-seaseq/master/wdl/util.wdl"
import "https://raw.githubusercontent.com/madetunj/wdl-seaseq/master/wdl/visualization.wdl" as viz
import "https://raw.githubusercontent.com/madetunj/wdl-seaseq/master/wdl/runspp.wdl"
import "https://raw.githubusercontent.com/madetunj/wdl-seaseq/master/wdl/sortbed.wdl"

workflow seaseq {
    input {
        File fastqfile
        File reference
        File reference_index
        File blacklistfile
        File chromsizes
        File gtffile
        Array[File]+ index_files
        Array[File]+ motif_databases
    }

    call fastqc.fastqc {
        input :
            inputfile=fastqfile
    }
    
    call util.basicfastqstats as bfs {
        input :
            fastqfile=fastqfile
    }
    
    call bowtie.bowtie {
        input :
            fastqfile=fastqfile,
            index_files=index_files,
            metricsfile=bfs.metrics_out
    }
    
    call samtools.viewsort {
        input :
            samfile=bowtie.samfile
    }
    
    call fastqc.fastqc as bamfqc {
        input :
            inputfile=viewsort.sortedbam
    }
    
    call samtools.indexstats {
        input :
            bamfile=viewsort.sortedbam
    }
    
    call bedtools.intersect as blacklist {
        input :
            fileA=viewsort.sortedbam,
            fileB=blacklistfile,
            default_location="BAM_files",
            nooverlap=true
    }
    
    call samtools.markdup {
        input :
            bamfile=blacklist.intersect_out
    }
    
    call samtools.indexstats as bklist {
        input :
            bamfile=blacklist.intersect_out
    }
    
    call samtools.indexstats as mkdup {
        input :
            bamfile=markdup.mkdupbam
    }
    
    call macs.macs {
        input :
            bamfile=blacklist.intersect_out
    }
    
    call macs.macs as all {
        input :
            bamfile=blacklist.intersect_out,
            keep_dup="all"
    }
    
    call macs.macs as nomodel {
        input :
            bamfile=blacklist.intersect_out,
            nomodel=true
    }
    
    call bamtogff.bamtogff {
        input :
            gtffile=gtffile,
            chromsizes=chromsizes,
            bamfile=markdup.mkdupbam,
            bamindex=mkdup.indexbam
    }
    
    call bedtools.bamtobed {
        input :
            bamfile=markdup.mkdupbam
    }
    
    call sicer.sicer {
        input :
            bedfile=bamtobed.bedfile
    }
    
    call motifs.motifs {
        input:
            reference=reference,
            reference_index=reference_index,
            bedfile=macs.peakbedfile,
            motif_databases=motif_databases
    }

    call util.flankbed {
        input :
            bedfile=macs.summitsfile,
            default_location="MOTIF_files"
    }
    
    call motifs.motifs as flank {
        input:
            reference=reference,
            reference_index=reference_index,
            bedfile=flankbed.flankbedfile,
            motif_databases=motif_databases
    }

    call rose.rose {
        input :
            gtffile=gtffile,
            bamfile=blacklist.intersect_out,
            bamindex=bklist.indexbam,
            bedfile_auto=macs.peakbedfile,
            bedfile_all=all.peakbedfile
    }

    call viz.visualization {
        input:
            wigfile=macs.wigfile,
            chromsizes=chromsizes,
            xlsfile=macs.peakxlsfile
    }
    
    call viz.visualization as vizall {
        input:
            wigfile=all.wigfile,
            chromsizes=chromsizes,
            xlsfile=all.peakxlsfile
    }
    
    call viz.visualization as viznomodel {
        input:
            wigfile=nomodel.wigfile,
            chromsizes=chromsizes,
            xlsfile=nomodel.peakxlsfile
    }

    call bedtools.bamtobed as tobed {
        input :
            bamfile=blacklist.intersect_out
    }
    
    call runspp.runspp {
        input:
            bamfile=blacklist.intersect_out
    }
    
    call sortbed.sortbed {
        input:
            bedfile=tobed.bedfile
    }
    
    call bedtools.intersect {
        input:
            fileA=macs.peakbedfile,
            fileB=sortbed.sortbed_out,
            countoverlap=true,
            sorted=true
    }
    
    call util.summarystats {
        input:
            bambed=tobed.bedfile,
            sppfile=runspp.spp_out,
            countsfile=intersect.intersect_out,
            peaksxls=macs.peakxlsfile,
            bamflag=indexstats.flagstats,
            rmdupflag=mkdup.flagstats,
            bkflag=bklist.flagstats,
            fastqczip=fastqc.zipfile,
            fastqmetrics=bfs.metrics_out,
            enhancers=rose.enhancers,
            superenhancers=rose.super_enhancers
    }
    
    output {
        #FASTQC
        File htmlfile = fastqc.htmlfile
        File zipfile = fastqc.zipfile
        File bam_htmlfile = bamfqc.htmlfile
        File bam_zipfile = bamfqc.zipfile

        #BASICMETRICS
        File metrics_out = bfs.metrics_out

        #FLANKBED
        File flankbedfile = flankbed.flankbedfile

        #BAMFILES
        File sortedbam = viewsort.sortedbam
        File mkdupbam = markdup.mkdupbam
        File bklistbam = blacklist.intersect_out
        File indexbam = indexstats.indexbam
        File bklist_indexbam = bklist.indexbam
        File mkdup_indexbam = mkdup.indexbam

        #MACS
        File peakbedfile = macs.peakbedfile
        File peakxlsfile = macs.peakxlsfile
        File summitsfile = macs.summitsfile
        File wigfile = macs.wigfile
        File all_peakbedfile = all.peakbedfile
        File all_peakxlsfile = all.peakxlsfile
        File all_summitsfile = all.summitsfile
        File all_wigfile = all.wigfile
        File nm_peakbedfile = nomodel.peakbedfile
        File nm_peakxlsfile = nomodel.peakxlsfile
        File nm_summitsfile = nomodel.summitsfile
        File nm_wigfile = nomodel.wigfile

        #SICER
        File scoreisland = sicer.scoreisland
        File sicer_wigfile = sicer.wigfile

        #ROSE
        File pngfile = rose.pngfile
        File? mapped_union = rose.mapped_union
        File? mapped_stitch = rose.mapped_stitch
        File enhancers = rose.enhancers
        File super_enhancers = rose.super_enhancers
        File? gff_file = rose.gff_file
        File? gff_union = rose.gff_union
        File? union_enhancers = rose.union_enhancers
        File? stitch_enhancers = rose.stitch_enhancers
        File? e_to_g_enhancers = rose.e_to_g_enhancers
        File? g_to_e_enhancers = rose.g_to_e_enhancers
        File? e_to_g_super_enhancers = rose.e_to_g_super_enhancers
        File? g_to_e_super_enhancers = rose.g_to_e_super_enhancers

        #MOTIFS
        File? ame_tsv = motifs.ame_tsv
        File? ame_html = motifs.ame_html
        File? ame_seq = motifs.ame_seq
        File meme = motifs.meme_out
        File meme_summary = motifs.meme_summary

        File? summit_ame_tsv = flank.ame_tsv
        File? summit_ame_html = flank.ame_html
        File? summit_ame_seq = flank.ame_seq
        File summit_meme = flank.meme_out
        File summit_meme_summary = flank.meme_summary

        #BAM2GFF
        File m_downstream = bamtogff.m_downstream
        File m_upstream = bamtogff.m_upstream
        File m_genebody = bamtogff.m_genebody
        File m_promoters = bamtogff.m_promoters
        File? pdf_gene = bamtogff.pdf_gene
        File? pdf_h_gene = bamtogff.pdf_h_gene
        File? png_h_gene = bamtogff.png_h_gene
        File? pdf_promoters = bamtogff.pdf_promoters
        File? pdf_h_promoters = bamtogff.pdf_h_promoters
        File? png_h_promoters = bamtogff.png_h_promoters

        #VISUALIZATION
        File bigwig = visualization.bigwig
        File norm_wig = visualization.norm_wig
        File tdffile = visualization.tdffile
        File n_bigwig = viznomodel.bigwig
        File n_norm_wig = viznomodel.norm_wig
        File n_tdffile = viznomodel.tdffile
        File a_bigwig = vizall.bigwig
        File a_norm_wig = vizall.norm_wig
        File a_tdffile = vizall.tdffile

        #QC-STATS
        File qc_statsfile = summarystats.statsfile
        File qc_htmlfile = summarystats.htmlfile
        File qc_textfile = summarystats.textfile
    }

}
