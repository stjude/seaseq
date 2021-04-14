version 1.0
import "https://raw.githubusercontent.com/stjude/seaseq/master/workflows/tasks/fastqc.wdl"
#import "https://raw.githubusercontent.com/stjude/seaseq/master/workflows/tasks/bedtools.wdl"
import "/Users/madetunj/Desktop/seaseq/workflows/tasks/bedtools.wdl"
#import "https://raw.githubusercontent.com/stjude/seaseq/master/workflows/tasks/bowtie.wdl"
import "/Users/madetunj/Desktop/seaseq/workflows/tasks/bowtie.wdl"
#import "https://raw.githubusercontent.com/stjude/seaseq/master/workflows/tasks/samtools.wdl"
import "/Users/madetunj/Desktop/seaseq/workflows/tasks/samtools.wdl"
import "https://raw.githubusercontent.com/stjude/seaseq/master/workflows/tasks/macs.wdl"
import "https://raw.githubusercontent.com/stjude/seaseq/master/workflows/tasks/bamtogff.wdl"
import "https://raw.githubusercontent.com/stjude/seaseq/master/workflows/tasks/sicer.wdl"
#import "https://raw.githubusercontent.com/stjude/seaseq/master/workflows/workflows/motifs.wdl"
import "/Users/madetunj/Desktop/seaseq/workflows/workflows/motifs.wdl"
import "https://raw.githubusercontent.com/stjude/seaseq/master/workflows/tasks/rose.wdl"
import "https://raw.githubusercontent.com/stjude/seaseq/master/workflows/tasks/util.wdl"
import "https://raw.githubusercontent.com/stjude/seaseq/master/workflows/workflows/visualization.wdl" as viz
import "https://raw.githubusercontent.com/stjude/seaseq/master/workflows/tasks/runspp.wdl"
import "https://raw.githubusercontent.com/stjude/seaseq/master/workflows/tasks/sortbed.wdl"
import "https://raw.githubusercontent.com/stjude/seaseq/master/workflows/tasks/sratoolkit.wdl" as sra

workflow seaseq {
    String pipeline_ver = 'v1.0.0'

    meta {
        title: 'SEAseq Analysis'
        summary: 'Single-End Antibody Sequencing (SEAseq) Pipeline'
        description: 'A comprehensive automated computational pipeline for all ChIP-Seq/CUT&RUN data analysis.'
        version: '1.0.0'
        details: {
            citation: 'pending',
            contactEmail: 'modupeore.adetunji@stjude.org',
            contactOrg: "St Jude Children's Research Hospital",
            contactUrl: "",
            upstreamLicenses: "MIT",
            upstreamUrl: 'https://github.com/stjude/seaseq',
            whatsNew: [
                {
                    version: "1.0",
                    changes: ["Initial release"]
                }
            ]
        }
        parameter_group: {
            reference_genome: {
                title: 'Reference genome',
                description: 'Genome specific files. e.g. reference FASTA, GTF, blacklist, motif databases, FASTA index, bowtie index .',
                help: 'Input reference genome files as defined. If some genome data are missing then analyses using such data will be skipped.'
            },
            input_genomic_data: {
                title: 'Input FASTQ data',
                description: 'Genomic input files for experiment.',
                help: 'Input one or more sample data and/or SRA identifiers.'
            } 
        }
    }
    input {
        # group: reference_genome
        File reference
        File? blacklist
        File gtf
        Array[File]? bowtie_index
        Array[File]? motif_databases
        
        # group: input_genomic_data
        Array[String]? sra_id 
        Array[File]? fastqfile

    }

    parameter_meta {
        reference: {
            description: 'Reference FASTA file',
            group: 'reference_genome',
            patterns: ["*.fa", "*.fasta", "*.fa.gz", "*.fasta.gz"]
        }
        blacklist: {
            description: 'Blacklist file in BED format',
            group: 'reference_genome',
            help: 'If defined, blacklist regions listed are excluded after reference alignment.',
            patterns: ["*.bed", "*.bed.gz"]
        }
        gtf: {
            description: 'gene annotation file (.gtf)',
            group: 'reference_genome',
            help: 'Input gene annotation file from RefSeq or GENCODE (.gtf).',
            patterns: ["*.gtf"]
        }
        bowtie_index: {
            description: 'bowtie v1 index files (*.ebwt)',
            group: 'reference_genome',
            help: 'If not defined, bowtie v1 index files are generated, will take a longer compute time.',
            patterns: ["*.ebwt"]
        }
        motif_databases: {
            description: 'One or more of the MEME suite motif databases (*.meme)',
            group: 'reference_genome',
            help: 'Input one or more motif databases available from the MEME suite (https://meme-suite.org/meme/db/motifs).',
            patterns: ["*.meme"]
        }
        sra_id: {
            description: 'One or more SRA (Sequence Read Archive) run identifiers',
            group: 'input_genomic_data',
            help: 'Input publicly available FASTQs (SRRs). Multiple SRRs are separated by commas (,).',
            example: 'SRR12345678'
        }
        fastqfile: {
            description: 'One or more FASTQs',
            group: 'input_genomic_data',
            help: 'Upload zipped FASTQ files.',
            patterns: ["*.fq.gz", "*.fastq.gz"]
        }
    }

    if ( defined(sra_id) ) {
        # download sample file(s) from SRA database
        # outputs:
        #    fastqdump.fastqfile : downloaded sample files in fastq.gz format 
        Array[String] fake_sra = [1] #buffer to allow for sra_id optionality
        Array[String] sra_id_ = select_first([sra_id, fake_sra])
        scatter (eachsra in sra_id_) {
            call sra.fastqdump {
                input :
                    sra_id=eachsra,
                    cloud="false"
            }
        }   

        Array[File] fastqfile_ = flatten(fastqdump.fastqfile)
    }

    if ( !defined(bowtie_index) ) {
        # create bowtie index when not provided
        call bowtie.index as bowtie_idx {
            input :
                reference=reference
        }
    }

    call samtools.faidx as samtools_faidx {
        input :
            reference=reference
    }

    # Merge files for the pipeline 
    Array[File] bowtie_index_ = flatten(select_all([bowtie_index, bowtie_idx.bowtie_indexes]))
    Array[File] fastqfiles = flatten(select_all([fastqfile, fastqfile_]))

    scatter (eachfastq in fastqfiles) {
        call fastqc.fastqc {
            input :
                inputfile=eachfastq,
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/QC/FastQC'
        }
        
        call util.basicfastqstats as bfs {
            input :
                fastqfile=eachfastq,
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/QC/SummaryStats'
        }
    
        call bowtie.bowtie {
            input :
                fastqfile=eachfastq,
                index_files=bowtie_index_,
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
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/QC/FastQC'
        }
    
        call samtools.indexstats {
            input :
                bamfile=viewsort.sortedbam,
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/BAM_files'
        }
    
        if (defined(blacklist)) {
            # remove blacklist regions
            String fake_blacklist = "" #buffer to allow for blacklist optionality
            File blacklist_ = select_first([blacklist, fake_blacklist])
            call bedtools.intersect as rmblklist {
                input :
                    fileA=viewsort.sortedbam,
                    fileB=blacklist_,
                    default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/BAM_files',
                    nooverlap=true
            }
            call samtools.indexstats as bklist {
                input :
                    bamfile=rmblklist.intersect_out,
                    default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/BAM_files'
            }
        }

        File downstream_bam = select_first([rmblklist.intersect_out, viewsort.sortedbam])
    
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
                pvalue = "1e-9",
                keep_dup="auto",
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/PEAKS/NARROW_peaks'+'/'+basename(downstream_bam,'\.bam') +'-p9_kd-auto'
        }
    
        call macs.macs as all {
            input :
                bamfile=downstream_bam,
                pvalue = "1e-9",
                keep_dup="all",
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/PEAKS/NARROW_peaks'+'/'+basename(downstream_bam,'\.bam') +'-p9_kd-all'
        }
        
        call macs.macs as nomodel {
            input :
                bamfile=downstream_bam,
                nomodel=true,
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/PEAKS/NARROW_peaks'+'/'+basename(downstream_bam,'\.bam') +'-nm'
        }
        
        call bamtogff.bamtogff {
            input :
                gtffile=gtf,
                chromsizes=samtools_faidx.chromsizes,
                bamfile=markdup.mkdupbam,
                bamindex=mkdup.indexbam,
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/BAM_Density'
        }
        
        call bedtools.bamtobed {
            input :
                bamfile=markdup.mkdupbam
        }
        
        call sicer.sicer {
            input :
                bedfile=bamtobed.bedfile,
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/PEAKS/BROAD_peaks'
        }
        
        if ( defined(motif_databases) ) {
            # motif prediction and enrichment analysis

            Array[String] fake_motif_databases = [1]
            Array[File] motif_databases_ = select_first([motif_databases, fake_motif_databases])
            call motifs.motifs {
                input:
                    reference=reference,
                    reference_index=samtools_faidx.faidx_file,
                    bedfile=macs.peakbedfile,
                    motif_databases=motif_databases_,
                    default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/MOTIFS'
            }
        
            call util.flankbed {
                input :
                    bedfile=macs.summitsfile,
                    default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/MOTIFS'
            }
            
            call motifs.motifs as flank {
                input:
                    reference=reference,
                    reference_index=samtools_faidx.faidx_file,
                    bedfile=flankbed.flankbedfile,
                    motif_databases=motif_databases_,
                    default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/MOTIFS'
            }
        }

        File rose_indexbam = select_first([bklist.indexbam, indexstats.indexbam])

        call rose.rose {
            input :
                gtffile=gtf,
                bamfile=downstream_bam,
                bamindex=rose_indexbam,
                bedfile_auto=macs.peakbedfile,
                bedfile_all=all.peakbedfile,
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/PEAKS/STITCHED_peaks'
        }

        call util.peaksanno {
            input :
                gtffile=gtf,
                bedfile=macs.peakbedfile,
                chromsizes=samtools_faidx.chromsizes,
                summitfile=macs.summitsfile,
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/PEAKS_Annotation'
        }

        call viz.visualization {
            input:
                wigfile=macs.wigfile,
                chromsizes=samtools_faidx.chromsizes,
                xlsfile=macs.peakxlsfile,
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/PEAKS_Display'
        }
        
        call viz.visualization as vizall {
            input:
                wigfile=all.wigfile,
                chromsizes=samtools_faidx.chromsizes,
                xlsfile=all.peakxlsfile,
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/PEAKS_Display'
        }
        
        call viz.visualization as viznomodel {
            input:
                wigfile=nomodel.wigfile,
                chromsizes=samtools_faidx.chromsizes,
                xlsfile=nomodel.peakxlsfile,
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/PEAKS_Display'
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
                default_location=sub(basename(eachfastq),'\.f.*q\.gz','')+'/QC/SummaryStats'
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

        #BAMFILES
        Array[File] sortedbam = viewsort.sortedbam
        Array[File] mkdupbam = markdup.mkdupbam
        Array[File?] bklistbam = rmblklist.intersect_out
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
        Array[File?] flankbedfile = flankbed.flankbedfile

        Array[File?] ame_tsv = motifs.ame_tsv
        Array[File?] ame_html = motifs.ame_html
        Array[File?] ame_seq = motifs.ame_seq
        Array[File?] meme = motifs.meme_out
        Array[File?] meme_summary = motifs.meme_summary

        Array[File?] summit_ame_tsv = flank.ame_tsv
        Array[File?] summit_ame_html = flank.ame_html
        Array[File?] summit_ame_seq = flank.ame_seq
        Array[File?] summit_meme = flank.meme_out
        Array[File?] summit_meme_summary = flank.meme_summary

        #BAM2GFF
        Array[File] m_downstream = bamtogff.m_downstream
        Array[File] m_upstream = bamtogff.m_upstream
        Array[File] m_genebody = bamtogff.m_genebody
        Array[File] m_promoters = bamtogff.m_promoters
        Array[File] densityplot = bamtogff.densityplot
        Array[File?] pdf_gene = bamtogff.pdf_gene
        Array[File?] pdf_h_gene = bamtogff.pdf_h_gene
        Array[File?] png_h_gene = bamtogff.png_h_gene
        Array[File?] pdf_promoters = bamtogff.pdf_promoters
        Array[File?] pdf_h_promoters = bamtogff.pdf_h_promoters
        Array[File?] png_h_promoters = bamtogff.png_h_promoters

        #PEAKS-ANNOTATION
        Array[File?] peak_promoters = peaksanno.peak_promoters
        Array[File?] peak_genebody = peaksanno.peak_genebody
        Array[File?] peak_window = peaksanno.peak_window
        Array[File?] peak_closest = peaksanno.peak_closest
        Array[File?] peak_comparison = peaksanno.peak_comparison
        Array[File?] gene_comparison = peaksanno.gene_comparison
        Array[File?] pdf_comparison = peaksanno.pdf_comparison

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
