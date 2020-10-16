version 1.0
import "https://raw.githubusercontent.com/stjude/seaseq/wdl-workflows/tasks/fastqc.wdl"
import "https://raw.githubusercontent.com/stjude/seaseq/wdl-workflows/tasks/bedtools.wdl"
import "https://raw.githubusercontent.com/stjude/seaseq/wdl-workflows/tasks/bowtie.wdl"
import "https://raw.githubusercontent.com/stjude/seaseq/wdl-workflows/tasks/samtools.wdl"
import "https://raw.githubusercontent.com/stjude/seaseq/wdl-workflows/tasks/macs.wdl"
import "https://raw.githubusercontent.com/stjude/seaseq/wdl-workflows/tasks/bamtogff.wdl"
import "https://raw.githubusercontent.com/stjude/seaseq/wdl-workflows/tasks/sicer.wdl"
import "https://raw.githubusercontent.com/stjude/seaseq/wdl-workflows/workflows/motifs.wdl"
import "https://raw.githubusercontent.com/stjude/seaseq/wdl-workflows/tasks/rose.wdl"
import "/home/madetunj/madetunj/seaseq/git-workflows/tasks/util.wdl"
#import "https://raw.githubusercontent.com/stjude/seaseq/wdl-workflows/tasks/util.wdl"
import "https://raw.githubusercontent.com/stjude/seaseq/wdl-workflows/workflows/visualization.wdl" as viz
import "https://raw.githubusercontent.com/stjude/seaseq/wdl-workflows/tasks/runspp.wdl"
import "https://raw.githubusercontent.com/stjude/seaseq/wdl-workflows/tasks/sortbed.wdl"
import "https://raw.githubusercontent.com/stjude/seaseq/wdl-workflows/tasks/sratoolkit.wdl" as sra

workflow seaseq {
    input {
        #input files 
        Array[String]? sra_id 
        Array[File]? fastqfile
        
        #reference genome + motif files
        File reference
        File? reference_index
        File? blacklistfile
        File chromsizes
        File gtffile
        Array[File]? index_files
        Array[File]+ motif_databases
        
        #additional variables
        String? gtf_feature = "gene"
    }

    if ( defined(sra_id) ) { 
        Array[String] fake_sra = [1] #buffer to allow for sra_id optionality
        Array[String] sra_id_ = select_first([sra_id, fake_sra])
        scatter (eachsra in sra_id_) {
            call sra.fastqdump {
                input :
                    sra_id=eachsra
            }
        }   

        Array[File] fastqfile_ = flatten(fastqdump.fastqfile)
    }

    if ( !defined(index_files) ) {
        call bowtie.index as bowtie_index {
            input:
                reference=reference
        }
    }

    if ( !defined(reference_index) ) {
        call samtools.faidx as samtools_faidx {
            input:
                reference=reference
        }
    }

    Array[File] index_files_ = flatten(select_all([index_files, bowtie_index.bowtie_indexes]))
    Array[File] fastqfiles = flatten(select_all([fastqfile, fastqfile_]))
    File reference_index_ = select_first([reference_index, samtools_faidx.faidx_file])

    scatter (eachfastq in fastqfiles) {
        call fastqc.fastqc {
            input :
                inputfile=eachfastq,
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/QC_files/FASTQC'
        }
        
        call util.basicfastqstats as bfs {
            input :
                fastqfile=eachfastq,
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/QC_files/STATS'
        }
    
        call bowtie.bowtie {
            input :
                fastqfile=eachfastq,
                index_files=index_files_,
                metricsfile=bfs.metrics_out
        }

        call samtools.viewsort {
            input :
                samfile=select_first(bowtie.samfile),
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/BAM_files'
        }
    
        call fastqc.fastqc as bamfqc {
            input :
                inputfile=viewsort.sortedbam,
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/QC_files/FASTQC'
        }
    
        call samtools.indexstats {
            input :
                bamfile=viewsort.sortedbam,
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/BAM_files'
        }
    
        if (defined(blacklistfile)) {
            String fake_blacklistfile = "" #buffer to allow for blacklistfile optionality
            File blacklistfile_ = select_first([blacklistfile, fake_blacklistfile])
            call bedtools.intersect as blacklist {
                input :
                    fileA=viewsort.sortedbam,
                    fileB=blacklistfile_,
                    default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/BAM_files',
                    nooverlap=true
            }
            call samtools.indexstats as bklist {
                input :
                    bamfile=blacklist.intersect_out,
                    default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/BAM_files'
            }
        }

        File downstream_bam = select_first([blacklist.intersect_out, viewsort.sortedbam])
    
        call samtools.markdup {
            input :
                bamfile=downstream_bam,
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/BAM_files'
        }

        call samtools.indexstats as mkdup {
            input :
                bamfile=markdup.mkdupbam,
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/BAM_files'
        }
    
        call macs.macs {
            input :
                bamfile=downstream_bam,
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/PEAKS_files/NARROW_peaks'
        }
    
        call macs.macs as all {
            input :
                bamfile=downstream_bam,
                keep_dup="all",
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/PEAKS_files/NARROW_peaks'
        }
        
        call macs.macs as nomodel {
            input :
                bamfile=downstream_bam,
                nomodel=true,
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/PEAKS_files/NARROW_peaks'
        }
        
        call bamtogff.bamtogff {
            input :
                feature=gtf_feature,
                gtffile=gtffile,
                chromsizes=chromsizes,
                bamfile=markdup.mkdupbam,
                bamindex=mkdup.indexbam,
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/BAMDensity_files'
        }
        
        call bedtools.bamtobed {
            input :
                bamfile=markdup.mkdupbam
        }
        
        call sicer.sicer {
            input :
                bedfile=bamtobed.bedfile,
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/PEAKS_files/BROAD_peaks'
        }
        
        call motifs.motifs {
            input:
                reference=reference,
                reference_index=reference_index_,
                bedfile=macs.peakbedfile,
                motif_databases=motif_databases,
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/MOTIF_files'
        }
    
        call util.flankbed {
            input :
                bedfile=macs.summitsfile,
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/MOTIF_files'
        }
        
        call motifs.motifs as flank {
            input:
                reference=reference,
                reference_index=reference_index_,
                bedfile=flankbed.flankbedfile,
                motif_databases=motif_databases,
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/MOTIF_files'
        }
        File rose_indexbam = select_first([bklist.indexbam, indexstats.indexbam])
        call rose.rose {
            input :
                feature=gtf_feature,
                gtffile=gtffile,
                bamfile=downstream_bam,
                bamindex=rose_indexbam,
                bedfile_auto=macs.peakbedfile,
                bedfile_all=all.peakbedfile,
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/PEAKS_files/STITCHED_REGIONS'
        }
    
        call viz.visualization {
            input:
                wigfile=macs.wigfile,
                chromsizes=chromsizes,
                xlsfile=macs.peakxlsfile,
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/PEAKDisplay_files'
        }
        
        call viz.visualization as vizall {
            input:
                wigfile=all.wigfile,
                chromsizes=chromsizes,
                xlsfile=all.peakxlsfile,
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/PEAKDisplay_files'
        }
        
        call viz.visualization as viznomodel {
            input:
                wigfile=nomodel.wigfile,
                chromsizes=chromsizes,
                xlsfile=nomodel.peakxlsfile,
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/PEAKDisplay_files'
        }
    
        call bedtools.bamtobed as tobed {
            input :
                bamfile=downstream_bam
        }
        
        call runspp.runspp {
            input:
                bamfile=downstream_bam
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
                superenhancers=rose.super_enhancers,
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/QC_files/STATS'
        }
    }
        
    output {
        #FASTQC
        Array[File] htmlfile = fastqc.htmlfile
        Array[File] zipfile = fastqc.zipfile
        Array[File] bam_htmlfile = bamfqc.htmlfile
        Array[File] bam_zipfile = bamfqc.zipfile

        #BASICMETRICS
        Array[File] metrics_out = bfs.metrics_out

        #FLANKBED
        Array[File] flankbedfile = flankbed.flankbedfile

        #BAMFILES
        Array[File] sortedbam = viewsort.sortedbam
        Array[File] mkdupbam = markdup.mkdupbam
        Array[File?] bklistbam = blacklist.intersect_out
        Array[File] indexbam = indexstats.indexbam
        Array[File?] bklist_indexbam = bklist.indexbam
        Array[File] mkdup_indexbam = mkdup.indexbam

        #MACS
        Array[File] peakbedfile = macs.peakbedfile
        Array[File] peakxlsfile = macs.peakxlsfile
        Array[File] summitsfile = macs.summitsfile
        Array[File] wigfile = macs.wigfile
        Array[File] all_peakbedfile = all.peakbedfile
        Array[File] all_peakxlsfile = all.peakxlsfile
        Array[File] all_summitsfile = all.summitsfile
        Array[File] all_wigfile = all.wigfile
        Array[File] nm_peakbedfile = nomodel.peakbedfile
        Array[File] nm_peakxlsfile = nomodel.peakxlsfile
        Array[File] nm_summitsfile = nomodel.summitsfile
        Array[File] nm_wigfile = nomodel.wigfile

        #SICER
        Array[File] scoreisland = sicer.scoreisland
        Array[File] sicer_wigfile = sicer.wigfile

        #ROSE
        Array[File] pngfile = rose.pngfile
        Array[File?] mapped_union = rose.mapped_union
        Array[File?] mapped_stitch = rose.mapped_stitch
        Array[File] enhancers = rose.enhancers
        Array[File] super_enhancers = rose.super_enhancers
        Array[File?] gff_file = rose.gff_file
        Array[File?] gff_union = rose.gff_union
        Array[File?] union_enhancers = rose.union_enhancers
        Array[File?] stitch_enhancers = rose.stitch_enhancers
        Array[File?] e_to_g_enhancers = rose.e_to_g_enhancers
        Array[File?] g_to_e_enhancers = rose.g_to_e_enhancers
        Array[File?] e_to_g_super_enhancers = rose.e_to_g_super_enhancers
        Array[File?] g_to_e_super_enhancers = rose.g_to_e_super_enhancers

        #MOTIFS
        Array[File?] ame_tsv = motifs.ame_tsv
        Array[File?] ame_html = motifs.ame_html
        Array[File?] ame_seq = motifs.ame_seq
        Array[File] meme = motifs.meme_out
        Array[File] meme_summary = motifs.meme_summary

        Array[File?] summit_ame_tsv = flank.ame_tsv
        Array[File?] summit_ame_html = flank.ame_html
        Array[File?] summit_ame_seq = flank.ame_seq
        Array[File] summit_meme = flank.meme_out
        Array[File] summit_meme_summary = flank.meme_summary

        #BAM2GFF
        Array[File] m_downstream = bamtogff.m_downstream
        Array[File] m_upstream = bamtogff.m_upstream
        Array[File] m_genebody = bamtogff.m_genebody
        Array[File] m_promoters = bamtogff.m_promoters
        Array[File?] pdf_gene = bamtogff.pdf_gene
        Array[File?] pdf_h_gene = bamtogff.pdf_h_gene
        Array[File?] png_h_gene = bamtogff.png_h_gene
        Array[File?] pdf_promoters = bamtogff.pdf_promoters
        Array[File?] pdf_h_promoters = bamtogff.pdf_h_promoters
        Array[File?] png_h_promoters = bamtogff.png_h_promoters

        #VISUALIZATION
        Array[File] bigwig = visualization.bigwig
        Array[File] norm_wig = visualization.norm_wig
        Array[File] tdffile = visualization.tdffile
        Array[File] n_bigwig = viznomodel.bigwig
        Array[File] n_norm_wig = viznomodel.norm_wig
        Array[File] n_tdffile = viznomodel.tdffile
        Array[File] a_bigwig = vizall.bigwig
        Array[File] a_norm_wig = vizall.norm_wig
        Array[File] a_tdffile = vizall.tdffile

        #QC-STATS
        Array[File] qc_statsfile = summarystats.statsfile
        Array[File] qc_htmlfile = summarystats.htmlfile
        Array[File] qc_textfile = summarystats.textfile
    }

}
