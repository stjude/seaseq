version 1.0
import "../tasks/fastqc.wdl"
import "../tasks/bedtools.wdl"
import "../tasks/bowtie.wdl"
import "../tasks/samtools.wdl"
import "../tasks/macs.wdl"
import "bamtogff.wdl"
import "../tasks/sicer.wdl"
import "motifs.wdl"
import "../tasks/rose.wdl"
import "../tasks/util.wdl"
import "visualization.wdl" as viz
import "../tasks/runspp.wdl"
import "../tasks/sortbed.wdl"

workflow paired_sample_analysis {
    String pipeline_ver = 'v2.0.0'

    meta {
        title: 'SEAseq Analysis'
        summary: 'Single-End Antibody Sequencing (SEAseq) Pipeline'
        description: 'A comprehensive automated computational pipeline for all ChIP-Seq/CUT&RUN data analysis.'
        version: '2.0.0'
        details: {
            citation: 'pending',
            contactEmail: 'modupeore.adetunji@stjude.org',
            contactOrg: "St Jude Children's Research Hospital",
            contactUrl: "",
            upstreamLicenses: "MIT",
            upstreamUrl: 'https://github.com/stjude/seaseq',
            whatsNew: [
                {
                    version: "2.0",
                    changes: ["version of case/sample + control", "single-end sequencing with input/control sequencing data", "Initial release"]
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
            },
            analysis_parameter: {
                title: 'Analysis parameter',
                description: 'Analysis settings needed for experiment.',
                help: 'Analysis settings; such output analysis file name.'
            }
        }
    }
    input {
        # group: reference_genome
        File reference
        File? blacklist
        File gtf
        Array[File]? motif_databases
        File chromsizes
        File faidx

        # group: input_genomic_data
        File sample_bam
        File sample_bai
        File control_bam
        File control_bai

        # group: analysis_parameter
        String? results_name
        Boolean run_motifs=true

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
            patterns: ["*.gtf", "*.gtf.gz", "*.gff", "*.gff.gz", "*.gff3", "*.gff3.gz"]
        }
        motif_databases: {
            description: 'One or more of the MEME suite motif databases (*.meme)',
            group: 'reference_genome',
            help: 'Input one or more motif databases available from the MEME suite (https://meme-suite.org/meme/db/motifs).',
            patterns: ["*.meme"]
        }
        sample_bam: {
            description: 'Sample BAM to analyze',
            group: 'input_genomic_data',
            help: 'Upload BAM file',
            patterns: ["*.bam"]
        }
        control_bam: {
            description: 'One input/control BAM',
            group: 'input_genomic_data',
            help: 'Upload BAM file.',
            patterns: ["*.bam"]
        }
        results_name: {
            description: 'Experiment results custom name',
            group: 'analysis_parameter',
            help: 'Input preferred analysis results name.',
            example: 'AllMerge_mapped'
        }
    }
    call macs.macs {
        input :
            bamfile=sample_bam,
            control=control_bam,
            pvalue = "1e-9",
            keep_dup="auto",
            output_name = if defined(results_name) then results_name + '-p9_kd-auto' else basename(sample_bam,'.bam') + '+control-p9_kd-auto',
            default_location = if defined(results_name) then results_name + '/PEAKS/NARROW_peaks/' + results_name + '-p9_kd-auto' else sub(basename(sample_bam),'.sorted.b.*$','') + '+control/PEAKS/NARROW_peaks/' + basename(sample_bam,'.bam') + '+control-p9_kd-auto'
    }

    call util.addreadme {
        input :
            default_location = if defined(results_name) then results_name + '/PEAKS' else sub(basename(sample_bam),'.sorted.b.*$','') + '+control/PEAKS'
    }
    
    call macs.macs as all {
        input :
            bamfile=sample_bam,
            control=control_bam,
            pvalue = "1e-9",
            keep_dup="all",
            output_name = if defined(results_name) then results_name + '-p9_kd-all' else basename(sample_bam,'.bam') + '+control-p9_kd-all',
            default_location = if defined(results_name) then results_name + '/PEAKS/NARROW_peaks/' + results_name + '-p9_kd-all' else sub(basename(sample_bam),'.sorted.b.*$','') + '+control/PEAKS/NARROW_peaks/' + basename(sample_bam,'.bam') + '+control-p9_kd-all'
    }

    call macs.macs as nomodel {
        input :
            bamfile=sample_bam,
            control=control_bam,
            nomodel=true,
            output_name = if defined(results_name) then results_name + '-nm' else basename(sample_bam,'.bam') + '+control-nm',
            default_location = if defined(results_name) then results_name + '/PEAKS/NARROW_peaks/' + results_name + '-nm' else sub(basename(sample_bam),'.sorted.b.*$','') + '+control/PEAKS/NARROW_peaks/' + basename(sample_bam,'.bam') + '+control-nm'
    }

    call bamtogff.bamtogff {
        input :
            gtffile=gtf,
            chromsizes=chromsizes,
            bamfile=sample_bam,
            bamindex=sample_bai,
            control_bamfile=control_bam,
            control_bamindex=control_bai,
            samplename=if defined(results_name) then results_name else basename(sample_bam,'.bam') + '+control',
            default_location=if defined(results_name) then results_name + '/BAM_Density' else sub(basename(sample_bam),'.sorted.b.*$','') + '+control/BAM_Density'
    }

    call bedtools.bamtobed as s_forsicerbed {
        input :
            bamfile=sample_bam
    }

    call bedtools.bamtobed as c_forsicerbed {
        input :
            bamfile=control_bam
    }

    call sicer.sicer {
        input :
            bedfile=s_forsicerbed.bedfile,
            control_bed=c_forsicerbed.bedfile,
            chromsizes=chromsizes,
            outputname=if defined(results_name) then results_name else basename(s_forsicerbed.bedfile,'.bed') + '+control',
            default_location=if defined(results_name) then results_name + '/PEAKS/BROAD_peaks' else sub(basename(sample_bam),'.sorted.b.*$','') + '+control/PEAKS/BROAD_peaks'
    }

    call rose.rose {
        input :
            gtffile=gtf,
            bamfile=sample_bam,
            bamindex=sample_bai,
            control=control_bam,
            controlindex=control_bai,
            bedfile_auto=macs.peakbedfile,
            bedfile_all=all.peakbedfile,
            default_location=if defined(results_name) then results_name + '/PEAKS/STITCHED_peaks' else sub(basename(sample_bam),'.sorted.b.*$','') + '+control/PEAKS/STITCHED_peaks'
    }

    call runspp.runspp {
        input:
            bamfile=sample_bam,
            control=control_bam
    }

    call runspp.runspp as only_s_runspp {
        input:
            bamfile=sample_bam
    }

    call runspp.runspp as only_c_runspp {
        input:
            bamfile=control_bam
    }

    String string_ctrlwig = ""
    call viz.visualization as c_visualization {
        input:
            wigfile=select_first([macs.ctrlwigfile, string_ctrlwig]),
            chromsizes=chromsizes,
            control=true,
            xlsfile=macs.peakxlsfile,
            default_location=if defined(results_name) then results_name + '/COVERAGE_files/NARROW_peaks/' + sub(basename(macs.peakbedfile),'_peaks.bed','') + '/control' else sub(basename(sample_bam),'.sorted.b.*$','') + '+control/COVERAGE_files/NARROW_peaks/' + sub(basename(macs.peakbedfile),'_peaks.bed','') + '/control'
    }

    call viz.visualization as c_vizall {
        input:
            wigfile=select_first([all.ctrlwigfile, string_ctrlwig]),
            chromsizes=chromsizes,
            control=true,
            xlsfile=all.peakxlsfile,
            default_location=if defined(results_name) then results_name + '/COVERAGE_files/NARROW_peaks/' + sub(basename(all.peakbedfile),'_peaks.bed','') + '/control' else sub(basename(sample_bam),'.sorted.b.*$','') + '+control/COVERAGE_files/NARROW_peaks/' + sub(basename(all.peakbedfile),'_peaks.bed','') + '/control'
    }
    call viz.visualization as c_viznomodel {
        input:
            wigfile=select_first([nomodel.ctrlwigfile, string_ctrlwig]),
            chromsizes=chromsizes,
            control=true,
            xlsfile=nomodel.peakxlsfile,
            default_location=if defined(results_name) then results_name + '/COVERAGE_files/NARROW_peaks/' + sub(basename(nomodel.peakbedfile),'_peaks.bed','') + '/control' else sub(basename(sample_bam),'.sorted.b.*$','') + '+control/COVERAGE_files/NARROW_peaks/' + sub(basename(nomodel.peakbedfile),'_peaks.bed','') + '/control'
    }

    call util.peaksanno {
        input :
            gtffile=gtf,
            bedfile=macs.peakbedfile,
            chromsizes=chromsizes,
            summitfile=macs.summitsfile,
            default_location=if defined(results_name) then results_name + '/PEAKS_Annotation/NARROW_peaks/' + sub(basename(macs.peakbedfile),'_peaks.bed','') else sub(basename(sample_bam),'.sorted.b.*$','') + '+control/PEAKS_Annotation/NARROW_peaks/' + sub(basename(macs.peakbedfile),'_peaks.bed','')
    }

    call util.peaksanno as all_peaksanno {
        input :
            gtffile=gtf,
            bedfile=all.peakbedfile,
            chromsizes=chromsizes,
            summitfile=all.summitsfile,
            default_location=if defined(results_name) then results_name + '/PEAKS_Annotation/NARROW_peaks/' + sub(basename(all.peakbedfile),'_peaks.bed','') else sub(basename(sample_bam),'.sorted.b.*$','') + '+control/PEAKS_Annotation/NARROW_peaks/' + sub(basename(all.peakbedfile),'_peaks.bed','')
    }

    call util.peaksanno as nomodel_peaksanno {
        input :
            gtffile=gtf,
            bedfile=nomodel.peakbedfile,
            chromsizes=chromsizes,
            summitfile=nomodel.summitsfile,
            default_location=if defined(results_name) then results_name + '/PEAKS_Annotation/NARROW_peaks/' + sub(basename(nomodel.peakbedfile),'_peaks.bed','') else sub(basename(sample_bam),'.sorted.b.*$','') + '+control/PEAKS_Annotation/NARROW_peaks/' + sub(basename(nomodel.peakbedfile),'_peaks.bed','')
    }

    call util.peaksanno as sicer_peaksanno {
        input :
            gtffile=gtf,
            bedfile=select_first([sicer.fdrisland, string_ctrlwig]),
            chromsizes=chromsizes,
            default_location=if defined(results_name) then results_name + '/PEAKS_Annotation/BROAD_peaks' else sub(basename(sample_bam),'.sorted.b.*$','') + '+control/PEAKS_Annotation/BROAD_peaks'
    }

    # Motif Analysis
    if (run_motifs) {
        call motifs.motifs {
            input:
                reference=reference,
                reference_index=faidx,
                bedfile=macs.peakbedfile,
                motif_databases=motif_databases,
                default_location=if defined(results_name) then results_name + '/MOTIFS' else sub(basename(sample_bam),'.sorted.b.*$','') + '+control/MOTIFS'
        }

        call util.flankbed {
            input :
                bedfile=macs.summitsfile,
                default_location=if defined(results_name) then results_name + '/MOTIFS' else sub(basename(sample_bam),'.sorted.b.*$','') + '+control/MOTIFS'
        }

        call motifs.motifs as flank {
            input:
                reference=reference,
                reference_index=faidx,
                bedfile=flankbed.flankbedfile,
                motif_databases=motif_databases,
                default_location=if defined(results_name) then results_name + '/MOTIFS' else sub(basename(sample_bam),'.sorted.b.*$','') + '+control/MOTIFS'
        }
    }

    call viz.visualization {
        input:
            wigfile=macs.wigfile,
            chromsizes=chromsizes,
            xlsfile=macs.peakxlsfile,
            default_location=if defined(results_name) then results_name + '/COVERAGE_files/NARROW_peaks/' + sub(basename(macs.peakbedfile),'_peaks.bed','') else sub(basename(sample_bam),'.sorted.b.*$','') + '+control/COVERAGE_files/NARROW_peaks/' + sub(basename(macs.peakbedfile),'_peaks.bed','')
    }

    call viz.visualization as vizall {
        input:
            wigfile=all.wigfile,
            chromsizes=chromsizes,
            xlsfile=all.peakxlsfile,
            default_location=if defined(results_name) then results_name + '/COVERAGE_files/NARROW_peaks/' + sub(basename(all.peakbedfile),'_peaks.bed','') else sub(basename(sample_bam),'.sorted.b.*$','') + '+control/COVERAGE_files/NARROW_peaks/' + sub(basename(all.peakbedfile),'_peaks.bed','')
    }

    call viz.visualization as viznomodel {
        input:
            wigfile=nomodel.wigfile,
            chromsizes=chromsizes,
            xlsfile=nomodel.peakxlsfile,
            default_location=if defined(results_name) then results_name + '/COVERAGE_files/NARROW_peaks/' + sub(basename(nomodel.peakbedfile),'_peaks.bed','') else sub(basename(sample_bam),'.sorted.b.*$','') + '+control/COVERAGE_files/NARROW_peaks/' + sub(basename(nomodel.peakbedfile),'_peaks.bed','')
    }

    call viz.visualization as vizsicer {
        input:
            wigfile=sicer.wigfile,
            chromsizes=chromsizes,
            default_location=if defined(results_name) then results_name + '/COVERAGE_files/BROAD_peaks' else sub(basename(sample_bam),'.sorted.b.*$','') + '+control/COVERAGE_files/BROAD_peaks'
    }

    #Peak Calling for Sample BAM only
    call macs.macs as only_s_macs {
        input :
            bamfile=sample_bam,
            pvalue = "1e-9",
            keep_dup="auto",
            default_location='SAMPLE/' + sub(basename(sample_bam),'.sorted.b.*$','') + '/PEAKS_forQC/' + basename(sample_bam,'.bam') + '-p9_kd-auto'
    }

    #Peak Calling for Control BAM only
    call macs.macs as only_c_macs {
        input :
            bamfile=control_bam,
            pvalue = "1e-9",
            keep_dup="auto",
            default_location='CONTROL/' + sub(basename(control_bam),'.sorted.b.*$','') + '/PEAKS_forQC/' + basename(control_bam,'.bam') + '-p9_kd-auto'
    }

    call bedtools.bamtobed as only_c_finalbed {
        input:
            bamfile=control_bam
    }

    call sortbed.sortbed as only_c_sortbed {
        input:
            bedfile=only_c_finalbed.bedfile
    }

    call bedtools.intersect as only_c_intersect {
        input:
            fileA=only_c_macs.peakbedfile,
            fileB=only_c_sortbed.sortbed_out,
            countoverlap=true,
            sorted=true
    }

    call bedtools.bamtobed as only_s_finalbed {
        input:
            bamfile=sample_bam
    }

    call sortbed.sortbed as only_s_sortbed {
        input:
            bedfile=only_s_finalbed.bedfile
    }

    call bedtools.intersect as only_s_intersect {
        input:
            fileA=only_s_macs.peakbedfile,
            fileB=only_s_sortbed.sortbed_out,
            countoverlap=true,
            sorted=true
    }

    call bedtools.intersect {
        input:
            fileA=macs.peakbedfile,
            fileB=only_s_sortbed.sortbed_out,
            countoverlap=true,
            sorted=true
    }

    output {
        #Bedtools
        File only_s_intersectfile = only_s_intersect.intersect_out
        File intersect_outfile = intersect.intersect_out
        File only_c_intersectfile = only_c_intersect.intersect_out
        File only_s_finalbedfile = only_s_finalbed.bedfile
        File only_c_finalbedfile = only_c_finalbed.bedfile

        #SPP
        File? only_s_runspp_file = only_s_runspp.spp_out
        File? only_c_runspp_file = only_c_runspp.spp_out
        File? runspp_file = runspp.spp_out

        #MACS
        File? peakbedfile = macs.peakbedfile
        File? peakxlsfile = macs.peakxlsfile
        File? summitsfile = macs.summitsfile
        File? negativexlsfile = macs.negativepeaks
        File? wigfile = macs.wigfile
        File? ctrlwigfile = macs.ctrlwigfile
        File? all_peakbedfile = all.peakbedfile
        File? all_peakxlsfile = all.peakxlsfile
        File? all_summitsfile = all.summitsfile
        File? all_negativexlsfile = all.negativepeaks
        File? all_wigfile = all.wigfile
        File? all_ctrlwigfile = all.ctrlwigfile
        File? nm_peakbedfile = nomodel.peakbedfile
        File? nm_peakxlsfile = nomodel.peakxlsfile
        File? nm_summitsfile = nomodel.summitsfile
        File? nm_negativexlsfile = nomodel.negativepeaks
        File? nm_wigfile = nomodel.wigfile
        File? nm_ctrlwigfile = nomodel.ctrlwigfile
        File? readme_peaks = addreadme.readme_peaks

        File? only_c_peakbedfile = only_c_macs.peakbedfile
        File? only_c_peakxlsfile = only_c_macs.peakxlsfile
        File? only_c_summitsfile = only_c_macs.summitsfile
        File? only_c_wigfile = only_c_macs.wigfile
        File? only_s_peakbedfile = only_s_macs.peakbedfile
        File? only_s_peakxlsfile = only_s_macs.peakxlsfile
        File? only_s_summitsfile = only_s_macs.summitsfile
        File? only_s_wigfile = only_s_macs.wigfile

        #SICER
        File? scoreisland = sicer.scoreisland
        File? sicer_wigfile = sicer.wigfile
        File? sicer_summary = sicer.summary
        File? sicer_fdrisland = sicer.fdrisland

        #ROSE
        File? pngfile = rose.pngfile
        File? mapped_union = rose.mapped_union
        File? mapped_stitch = rose.mapped_stitch
        File? enhancers = rose.enhancers
        File? super_enhancers = rose.super_enhancers
        File? gff_file = rose.gff_file
        File? gff_union = rose.gff_union
        File? union_enhancers = rose.union_enhancers
        File? stitch_enhancers = rose.stitch_enhancers
        File? e_to_g_enhancers = rose.e_to_g_enhancers
        File? g_to_e_enhancers = rose.g_to_e_enhancers
        File? e_to_g_super_enhancers = rose.e_to_g_super_enhancers
        File? g_to_e_super_enhancers = rose.g_to_e_super_enhancers

        #MOTIFS
        File? flankbedfile = flankbed.flankbedfile

        File? ame_tsv = motifs.ame_tsv
        File? ame_html = motifs.ame_html
        File? ame_seq = motifs.ame_seq
        File? meme = motifs.meme_out
        File? meme_summary = motifs.meme_summary

        File? summit_ame_tsv = flank.ame_tsv
        File? summit_ame_html = flank.ame_html
        File? summit_ame_seq = flank.ame_seq
        File? summit_meme = flank.meme_out
        File? summit_meme_summary = flank.meme_summary

        #BAM2GFF
        File? s_matrices = bamtogff.s_matrices
        File? c_matrices = bamtogff.c_matrices
        File? densityplot = bamtogff.densityplot
        File? pdf_gene = bamtogff.pdf_gene
        File? pdf_h_gene = bamtogff.pdf_h_gene
        File? png_h_gene = bamtogff.png_h_gene
        File? pdf_promoters = bamtogff.pdf_promoters
        File? pdf_h_promoters = bamtogff.pdf_h_promoters
        File? png_h_promoters = bamtogff.png_h_promoters

        #PEAKS-ANNOTATION
        File? peak_promoters = peaksanno.peak_promoters
        File? peak_genebody = peaksanno.peak_genebody
        File? peak_window = peaksanno.peak_window
        File? peak_closest = peaksanno.peak_closest
        File? peak_comparison = peaksanno.peak_comparison
        File? gene_comparison = peaksanno.gene_comparison
        File? pdf_comparison = peaksanno.pdf_comparison

        File? all_peak_promoters = all_peaksanno.peak_promoters
        File? all_peak_genebody = all_peaksanno.peak_genebody
        File? all_peak_window = all_peaksanno.peak_window
        File? all_peak_closest = all_peaksanno.peak_closest
        File? all_peak_comparison = all_peaksanno.peak_comparison
        File? all_gene_comparison = all_peaksanno.gene_comparison
        File? all_pdf_comparison = all_peaksanno.pdf_comparison

        File? nomodel_peak_promoters = nomodel_peaksanno.peak_promoters
        File? nomodel_peak_genebody = nomodel_peaksanno.peak_genebody
        File? nomodel_peak_window = nomodel_peaksanno.peak_window
        File? nomodel_peak_closest = nomodel_peaksanno.peak_closest
        File? nomodel_peak_comparison = nomodel_peaksanno.peak_comparison
        File? nomodel_gene_comparison = nomodel_peaksanno.gene_comparison
        File? nomodel_pdf_comparison = nomodel_peaksanno.pdf_comparison

        File? sicer_peak_promoters = sicer_peaksanno.peak_promoters
        File? sicer_peak_genebody = sicer_peaksanno.peak_genebody
        File? sicer_peak_window = sicer_peaksanno.peak_window
        File? sicer_peak_closest = sicer_peaksanno.peak_closest
        File? sicer_peak_comparison = sicer_peaksanno.peak_comparison
        File? sicer_gene_comparison = sicer_peaksanno.gene_comparison
        File? sicer_pdf_comparison = sicer_peaksanno.pdf_comparison

        #VISUALIZATION
        File? bigwig = visualization.bigwig
        File? norm_wig = visualization.norm_wig
        File? tdffile = visualization.tdffile
        File? n_bigwig = viznomodel.bigwig
        File? n_norm_wig = viznomodel.norm_wig
        File? n_tdffile = viznomodel.tdffile
        File? a_bigwig = vizall.bigwig
        File? a_norm_wig = vizall.norm_wig
        File? a_tdffile = vizall.tdffile

        File? c_bigwig = c_visualization.bigwig
        File? c_norm_wig = c_visualization.norm_wig
        File? c_tdffile = c_visualization.tdffile
        File? c_n_bigwig = c_viznomodel.bigwig
        File? c_n_norm_wig = c_viznomodel.norm_wig
        File? c_n_tdffile = c_viznomodel.tdffile
        File? c_a_bigwig = c_vizall.bigwig
        File? c_a_norm_wig = c_vizall.norm_wig
        File? c_a_tdffile = c_vizall.tdffile

        File? s_bigwig = vizsicer.bigwig
        File? s_norm_wig = vizsicer.norm_wig
        File? s_tdffile = vizsicer.tdffile
    }


}
