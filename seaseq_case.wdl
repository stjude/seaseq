version 1.0
import "workflows/tasks/fastqc.wdl"
import "workflows/tasks/bedtools.wdl"
import "workflows/tasks/bowtie.wdl"
import "workflows/tasks/samtools.wdl"
import "workflows/tasks/macs.wdl"
import "workflows/workflows/bamtogff.wdl"
import "workflows/tasks/sicer.wdl"
import "workflows/workflows/motifs.wdl"
import "workflows/tasks/rose.wdl"
import "workflows/tasks/util.wdl"
import "workflows/workflows/visualization.wdl" as viz
import "workflows/workflows/mapping.wdl"
import "workflows/tasks/runspp.wdl"
import "workflows/tasks/sortbed.wdl"
import "workflows/tasks/sratoolkit.wdl" as sra

workflow seaseq {
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
                    changes: ["version of case/sample only", "single-end sequencing with input/control sequencing data", "Initial release"]
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
        Array[File]? bowtie_index
        Array[File]? motif_databases

        # group: input_genomic_data
        Array[String]? sample_sraid
        Array[File]? sample_fastq

        # group: analysis_parameter
        String? sample_name

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
        sample_sraid: {
            description: 'One or more sample SRA (Sequence Read Archive) run identifiers',
            group: 'input_genomic_data',
            help: 'Input publicly available FASTQs (SRRs). Multiple SRRs are separated by commas (,).',
            example: 'SRR12345678'
        }
        sample_fastq: {
            description: 'One or more sample FASTQs',
            group: 'input_genomic_data',
            help: 'Upload zipped FASTQ files.',
            patterns: ["*fq", "*.fq.gz", "*.fastq", "*.fastq.gz"]
        }
        sample_name: {
            description: 'Sample results custom name',
            group: 'analysis_parameter',
            help: 'Input preferred analysis results file name (recommended if multiple FASTQs are provided).',
            example: 'AllMerge_mapped'
        }
    }

### ---------------------------------------- ###
### ------------ S E C T I O N 1 ----------- ###
### ------ pre-process analysis files ------ ###
### ---------------------------------------- ###

    if ( defined(sample_sraid) ) {
        # Download sample file(s) from SRA database
        # outputs:
        #    fastqdump.fastqfile : downloaded sample files in fastq.gz format
        Array[String] string_sra = [1] #buffer to allow for sra_id optionality
        Array[String] s_sraid = select_first([sample_sraid, string_sra])
        scatter (eachsra in s_sraid) {
            call sra.fastqdump {
                input :
                    sra_id=eachsra,
                    cloud="false"
            }
        }

        Array[File] sample_fastqfile_ = flatten(fastqdump.fastqfile)
    }

    #Generating INDEX files
    #1. Bowtie INDEX files if not provided
    if ( !defined(bowtie_index) ) {
        # create bowtie index when not provided
        call bowtie.index as bowtie_idx {
            input :
                reference=reference
        }
    }

    if ( defined(bowtie_index) ) {
        # check total number of bowtie indexes provided
        Array[String] string_bowtie_index = [1] #buffer to allow for bowtie_index optionality
        Array[File] int_bowtie_index = select_first([bowtie_index, string_bowtie_index])
        if ( length(int_bowtie_index) != 6 ) {
            # create bowtie index if 6 index files aren't provided
            call bowtie.index as bowtie_idx_2 {
                input :
                    reference=reference
            }
        }
    }

    call samtools.faidx as samtools_faidx {
        # create FASTA index and chrom sizes files
        input :
            reference=reference
    }


    Array[File] bowtie_index_ = select_first([bowtie_idx_2.bowtie_indexes, bowtie_idx.bowtie_indexes, bowtie_index])
    Array[File] s_fastqfiles = flatten(select_all([sample_fastq, sample_fastqfile_]))

### ------------------------------------------------- ###
### ---------------- S E C T I O N 2 ---------------- ###
### ---- A: analysis if multiple FASTQs provided ---- ###
### ------------------------------------------------- ###

    # if multiple fastqfiles are provided
    Boolean multi_sample_fastq = if length(s_fastqfiles) > 1 then true else false
    Boolean one_sample_fastq = if length(s_fastqfiles) == 1 then true else false

    if ( multi_sample_fastq ) {
        scatter (eachfastq in s_fastqfiles) {
            # Execute analysis on each fastq file provided
            # Analysis executed:
            #   FastQC
            #   FASTQ read length distribution
            #   Reference Alignment using Bowtie (-k2 -m2)
            #   Convert SAM to BAM
            #   FastQC on BAM files
            #   Remove Blacklists (if provided)
            #   Remove read duplicates
            #   Summary statistics on FASTQs
            #   Combine html files into one for easy viewing
            
            call fastqc.fastqc as fastqc {
                input :
                    inputfile=eachfastq,
                    default_location='SAMPLE/' + sub(basename(eachfastq),'\.f.*q\.gz','') + '/QC/FastQC'
            }

            call util.basicfastqstats as bfs {
                input :
                    fastqfile=eachfastq,
                    default_location='SAMPLE/' + sub(basename(eachfastq),'\.f.*q\.gz','') + '/QC/SummaryStats'
            }

            call mapping.mapping as each_mapping {
                input :
                    fastqfile=eachfastq,
                    index_files=bowtie_index_,
                    metricsfile=bfs.metrics_out,
                    blacklist=blacklist,
                    default_location='SAMPLE/' + sub(basename(eachfastq),'\.f.*q\.gz','') + '/BAM_files'
            }

            call fastqc.fastqc as bamfqc {
                input :
                    inputfile=each_mapping.sorted_bam,
                    default_location='SAMPLE/' + sub(basename(eachfastq),'\.f.*q\.gz','') + '/QC/FastQC'
            }

            call runspp.runspp as indv_runspp {
                input:
                    bamfile=select_first([each_mapping.bklist_bam, each_mapping.sorted_bam])
            }

            call bedtools.bamtobed as indv_bamtobed {
                input:
                    bamfile=select_first([each_mapping.bklist_bam, each_mapping.sorted_bam])
            }

            call util.evalstats as s_summarystats {
                input:
                    sppfile=indv_runspp.spp_out,
                    bambed=indv_bamtobed.bedfile,
                    sample_bamflag=each_mapping.bam_stats,
                    sample_rmdupflag=each_mapping.mkdup_stats,
                    sample_bkflag=each_mapping.bklist_stats,
                    sample_fastqczip=fastqc.zipfile,
                    fastqmetrics=bfs.metrics_out,
                    default_location='SAMPLE/' + sub(basename(eachfastq),'\.f.*q\.gz','') + '/QC/SummaryStats'
            }
        } # end scatter (for each sample fastq)

        # MERGE BAM FILES
        # Execute analysis on merge bam file
        # Analysis executed:
        #   Merge BAM (if more than 1 fastq is provided)
        #   FastQC on Merge BAM (AllMerge_<number>_mapped)

        # merge bam files and perform fasTQC if more than one is provided
        call util.mergehtml {
            input:
                htmlfiles=s_summarystats.htmlfile,
                txtfiles=s_summarystats.textfile,
                default_location='SAMPLE',
                outputfile = 'AllMapped_' + length(s_fastqfiles) + '_seaseq-summary-stats.html'
        }

        call samtools.mergebam {
            input:
                bamfiles=each_mapping.sorted_bam,
                default_location = if defined(sample_name) then sample_name + '/BAM_files' else 'AllMerge_' + length(each_mapping.sorted_bam) + '_mapped' + '/BAM_files',
                outputfile = if defined(sample_name) then sample_name + '.sorted.bam' else 'AllMerge_' + length(s_fastqfiles) + '_mapped.sorted.bam'
        }

        call fastqc.fastqc as mergebamfqc {
            input:
	        inputfile=mergebam.mergebam,
                default_location=sub(basename(mergebam.mergebam),'\.sorted\.b.*$','') + '/QC/FastQC'
        }

        call samtools.indexstats as mergeindexstats {
            input:
                bamfile=mergebam.mergebam,
                default_location=sub(basename(mergebam.mergebam),'\.sorted\.b.*$','') + '/BAM_files'
        }

        if ( defined(blacklist) ) {
            # remove blacklist regions
            String string_blacklist = "" #buffer to allow for blacklist optionality
            File blacklist_ = select_first([blacklist, string_blacklist])
            call bedtools.intersect as merge_rmblklist {
                input :
                    fileA=mergebam.mergebam,
                    fileB=blacklist_,
                    default_location=sub(basename(mergebam.mergebam),'\.sorted\.b.*$','') + '/BAM_files',
                    nooverlap=true
            }
            call samtools.indexstats as merge_bklist {
                input :
                    bamfile=merge_rmblklist.intersect_out,
                    default_location=sub(basename(mergebam.mergebam),'\.sorted\.b.*$','') + '/BAM_files'
            }
        } # end if blacklist provided

        File mergebam_afterbklist = select_first([merge_rmblklist.intersect_out, mergebam.mergebam])

        call samtools.markdup as merge_markdup {
            input :
                bamfile=mergebam_afterbklist,
                default_location=sub(basename(mergebam_afterbklist),'\.sorted\.b.*$','') + '/BAM_files'
        }

        call samtools.indexstats as merge_mkdup {
            input :
                bamfile=merge_markdup.mkdupbam,
                default_location=sub(basename(mergebam_afterbklist),'\.sorted\.b.*$','') + '/BAM_files'
        }
    } # end if length(fastqfiles) > 1: multi_sample_fastq

### ---------------------------------------- ###
### ------------ S E C T I O N 2 ----------- ###
### -- B: analysis if one FASTQ provided --- ###
### ---------------------------------------- ###

    # if only one fastqfile is provided
    if ( one_sample_fastq ) {
        # Execute analysis on each fastq file provided
        # Analysis executed:
        #   FastQC
        #   FASTQ read length distribution
        #   Reference Alignment using Bowtie (-k2 -m2)
        #   Convert SAM to BAM
        #   FastQC on BAM files
        #   Remove Blacklists (if provided)
        #   Remove read duplicates
        #   Summary statistics on FASTQs
        #   Combine html files into one for easy viewing
        if ( defined(sample_name) ) {
            String string_prefix = "" #buffer to allow for sample_name optionality
            String sample_name_ = select_first([sample_name, string_prefix])
            call util.linkname {
                input:
                    prefix=sample_name_,
                    in_fastq=s_fastqfiles[0]
            }
        }

        File uno_s_fastqfile = select_first([linkname.out_fastq, s_fastqfiles[0]])

        call fastqc.fastqc as uno_s_fastqc {
            input :
                inputfile=uno_s_fastqfile,
                default_location=sub(basename(uno_s_fastqfile),'\.f.*q\.gz','') + '/QC/FastQC'
        }

        call util.basicfastqstats as uno_s_bfs {
            input :
                fastqfile=uno_s_fastqfile,
                default_location=sub(basename(uno_s_fastqfile),'\.f.*q\.gz','') + '/QC/SummaryStats'
        }

        call mapping.mapping {
            input :
                fastqfile=uno_s_fastqfile,
                index_files=bowtie_index_,
                metricsfile=uno_s_bfs.metrics_out,
                blacklist=blacklist,
                default_location=sub(basename(uno_s_fastqfile),'\.f.*q\.gz','') + '/BAM_files'
        }

        call fastqc.fastqc as uno_s_bamfqc {
            input :
                inputfile=mapping.sorted_bam,
                default_location=sub(basename(uno_s_fastqfile),'\.f.*q\.gz','') + '/QC/FastQC'
        }

        call runspp.runspp as uno_s_runspp {
            input:
                bamfile=select_first([mapping.bklist_bam, mapping.sorted_bam])
        }

        call bedtools.bamtobed as uno_s_bamtobed {
            input:
                bamfile=select_first([mapping.bklist_bam, mapping.sorted_bam])
        }

        call util.evalstats as uno_s_summarystats {
            input:
                sppfile=uno_s_runspp.spp_out,
                bambed=uno_s_bamtobed.bedfile,
                sample_bamflag=mapping.bam_stats,
                sample_rmdupflag=mapping.mkdup_stats,
                sample_bkflag=mapping.bklist_stats,
                sample_fastqczip=uno_s_fastqc.zipfile,
                fastqmetrics=uno_s_bfs.metrics_out,
                outputfile = sub(basename(s_fastqfiles[0]),'\.f.*q\.gz', '-stats.csv'),
                outputhtml = sub(basename(s_fastqfiles[0]),'\.f.*q\.gz', '-stats.html'),
                outputtext = sub(basename(s_fastqfiles[0]),'\.f.*q\.gz', '-stats.txt'),
                configml = sub(basename(s_fastqfiles[0]),'\.f.*q\.gz', '-config.ml'),
                default_location='SAMPLE/' + sub(basename(s_fastqfiles[0]),'\.f.*q\.gz','') + '/QC/SummaryStats'
        }
    } # end if length(fastqfiles) == 1: one_sample_fastq

### ---------------------------------------- ###
### ------------ S E C T I O N 3 ----------- ###
### ----------- ChIP-seq analysis ---------- ###
### ---------------------------------------- ###

    # ChIP-seq and downstream analysis
    # Execute analysis on merge bam file
    # Analysis executed:
    #   FIRST: Check if reads are mapped
    #   Peaks identification (SICER, MACS, ROSE)
    #   Motif analysis
    #   Complete Summary statistics

    #collate correct files for downstream analysis
    File sample_bam = select_first([mergebam_afterbklist, mapping.bklist_bam, mapping.sorted_bam])

    call macs.macs {
        input :
            bamfile=sample_bam,
            pvalue = "1e-9",
            keep_dup="auto",
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/PEAKS/NARROW_peaks' + '/' + basename(sample_bam,'\.bam') + '-p9_kd-auto'
    }

    call macs.macs as all {
        input :
            bamfile=sample_bam,
            pvalue = "1e-9",
            keep_dup="all",
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/PEAKS/NARROW_peaks' + '/' + basename(sample_bam,'\.bam') + '-p9_kd-all'
    }

    call macs.macs as nomodel {
        input :
            bamfile=sample_bam,
            nomodel=true,
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/PEAKS/NARROW_peaks' + '/' + basename(sample_bam,'\.bam') + '-nm'
    }

    call bamtogff.bamtogff {
        input :
            gtffile=gtf,
            chromsizes=samtools_faidx.chromsizes,
            bamfile=select_first([merge_markdup.mkdupbam, mapping.mkdup_bam]),
            bamindex=select_first([merge_mkdup.indexbam, mapping.mkdup_index]),
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/BAM_Density'
    }

    call bedtools.bamtobed as forsicerbed {
        input :
            bamfile=select_first([merge_markdup.mkdupbam, mapping.mkdup_bam])
    }
    
    call sicer.sicer {
        input :
            bedfile=forsicerbed.bedfile,
            chromsizes=samtools_faidx.chromsizes,
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/PEAKS/BROAD_peaks'
    }

    call rose.rose {
        input :
            gtffile=gtf,
            bamfile=sample_bam,
            bamindex=select_first([merge_bklist.indexbam, mergeindexstats.indexbam, mapping.bklist_index, mapping.bam_index]),
            bedfile_auto=macs.peakbedfile,
            bedfile_all=all.peakbedfile,
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/PEAKS/STITCHED_peaks'
    }

    call runspp.runspp {
        input:
            bamfile=sample_bam
    }

    call util.peaksanno {
        input :
            gtffile=gtf,
            bedfile=macs.peakbedfile,
            chromsizes=samtools_faidx.chromsizes,
            summitfile=macs.summitsfile,
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/PEAKS_Annotation/NARROW_peaks' + '/' + sub(basename(macs.peakbedfile),'\_peaks.bed','')
    }

    call util.peaksanno as all_peaksanno {
        input :
            gtffile=gtf,
            bedfile=all.peakbedfile,
            chromsizes=samtools_faidx.chromsizes,
            summitfile=all.summitsfile,
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/PEAKS_Annotation/NARROW_peaks' + '/' + sub(basename(all.peakbedfile),'\_peaks.bed','')
    }

    call util.peaksanno as nomodel_peaksanno {
        input :
            gtffile=gtf,
            bedfile=nomodel.peakbedfile,
            chromsizes=samtools_faidx.chromsizes,
            summitfile=nomodel.summitsfile,
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/PEAKS_Annotation/NARROW_peaks' + '/' + sub(basename(nomodel.peakbedfile),'\_peaks.bed','')
    }

    call util.peaksanno as sicer_peaksanno {
        input :
            gtffile=gtf,
            bedfile=sicer.scoreisland,
            chromsizes=samtools_faidx.chromsizes,
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/PEAKS_Annotation/BROAD_peaks'
    }

    # Motif Analysis
    call motifs.motifs {
        input:
            reference=reference,
            reference_index=samtools_faidx.faidx_file,
            bedfile=macs.peakbedfile,
            motif_databases=motif_databases,
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/MOTIFS'
    }

    call util.flankbed {
        input :
            bedfile=macs.summitsfile,
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/MOTIFS'
    }

    call motifs.motifs as flank {
        input:
            reference=reference,
            reference_index=samtools_faidx.faidx_file,
            bedfile=flankbed.flankbedfile,
            motif_databases=motif_databases,
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/MOTIFS'
    }

    call viz.visualization {
        input:
            wigfile=macs.wigfile,
            chromsizes=samtools_faidx.chromsizes,
            xlsfile=macs.peakxlsfile,
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/COVERAGE_files/NARROW_peaks' + '/' + sub(basename(macs.peakbedfile),'\_peaks.bed','')
    }

    call viz.visualization as vizall {
        input:
            wigfile=all.wigfile,
            chromsizes=samtools_faidx.chromsizes,
            xlsfile=all.peakxlsfile,
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/COVERAGE_files/NARROW_peaks' + '/' + sub(basename(all.peakbedfile),'\_peaks.bed','')
    }

    call viz.visualization as viznomodel {
        input:
            wigfile=nomodel.wigfile,
            chromsizes=samtools_faidx.chromsizes,
            xlsfile=nomodel.peakxlsfile,
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/COVERAGE_files/NARROW_peaks' + '/' + sub(basename(nomodel.peakbedfile),'\_peaks.bed','')
    }

    call viz.visualization as vizsicer {
        input:
            wigfile=sicer.wigfile,
            chromsizes=samtools_faidx.chromsizes,
            default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/COVERAGE_files/BROAD_peaks'
    }

    call bedtools.bamtobed as finalbed {
        input:
            bamfile=sample_bam
    }

    call sortbed.sortbed {
        input:
            bedfile=finalbed.bedfile
    }

    call bedtools.intersect {
        input:
            fileA=macs.peakbedfile,
            fileB=sortbed.sortbed_out,
            countoverlap=true,
            sorted=true
    }

### ---------------------------------------- ###
### ------------ S E C T I O N 4 ----------- ###
### ---------- Summary Statistics ---------- ###
### ---------------------------------------- ###

    String string_qual = "" #buffer to allow for optionality in if statement

    #SUMMARY STATISTICS
    if ( one_sample_fastq ) {
        call util.evalstats as uno_summarystats {
            # SUMMARY STATISTICS of sample file (only 1 sample file provided)
            input:
                bambed=finalbed.bedfile,
                sppfile=runspp.spp_out,
                sample_fastqczip=select_first([uno_s_bamfqc.zipfile, string_qual]),
                sample_bamflag=mapping.bam_stats,
                sample_rmdupflag=mapping.mkdup_stats,
                sample_bkflag=mapping.bklist_stats,
                countsfile=intersect.intersect_out,
                peaksxls=macs.peakxlsfile,
                enhancers=rose.enhancers,
                superenhancers=rose.super_enhancers,
                default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/QC/SummaryStats'
        }
    } # end if one_fastq

    if ( multi_sample_fastq ) {
        call util.evalstats as merge_summarystats {
            # SUMMARY STATISTICS of all samples files (more than 1 sample file provided)
            input:
                bambed=finalbed.bedfile,
                sppfile=runspp.spp_out,
                sample_fastqczip=select_first([mergebamfqc.zipfile, string_qual]),
                sample_bamflag=mergeindexstats.flagstats,
                sample_rmdupflag=merge_mkdup.flagstats,
                sample_bkflag=merge_bklist.flagstats,
                countsfile=intersect.intersect_out,
                peaksxls=macs.peakxlsfile,
                enhancers=rose.enhancers,
                superenhancers=rose.super_enhancers,
                default_location=sub(basename(sample_bam),'\.sorted\.b.*$','') + '/QC/SummaryStats'
        }
    } # end if multi_fastq

    call util.summaryreport as overallsummary {
        # Presenting all quality stats for the analysis
        input:
            sampleqc_html=select_first([uno_s_summarystats.htmlfile, mergehtml.mergefile]),
            overallqc_html=select_first([uno_summarystats.htmlfile, merge_summarystats.htmlfile]),
            sampleqc_txt=select_first([uno_s_summarystats.textfile, mergehtml.mergetxt]),
            overallqc_txt=select_first([uno_summarystats.textfile, merge_summarystats.textfile])
    }

    output {
        #FASTQC
        Array[File?]? htmlfile = fastqc.htmlfile
        Array[File?]? zipfile = fastqc.zipfile
        Array[File?]? bam_htmlfile = bamfqc.htmlfile
        Array[File?]? bam_zipfile = bamfqc.zipfile
        File? mergebam_htmlfile = mergebamfqc.htmlfile
        File? mergebam_zipfile = mergebamfqc.zipfile
        File? uno_htmlfile = uno_s_fastqc.htmlfile
        File? uno_zipfile = uno_s_fastqc.zipfile
        File? uno_bam_htmlfile = uno_s_bamfqc.htmlfile
        File? uno_bam_zipfile = uno_s_bamfqc.zipfile

        #BASICMETRICS
        Array[File?]? metrics_out = bfs.metrics_out
        File? uno_metrics_out = uno_s_bfs.metrics_out

        #BAMFILES
        Array[File?]? sortedbam = each_mapping.sorted_bam
        Array[File?]? indexbam = each_mapping.bam_index
        Array[File?]? indv_bkbam = each_mapping.bklist_bam
        Array[File?]? indv_bkindexbam = each_mapping.bklist_index
        Array[File?]? indv_rmbam = each_mapping.mkdup_bam
        Array[File?]? indv_rmindexbam = each_mapping.mkdup_index

        File? uno_sortedbam = mapping.sorted_bam
        File? uno_indexstatsbam = mapping.bam_index
        File? uno_bkbam = mapping.bklist_bam
        File? uno_bkindexbam = mapping.bklist_index
        File? uno_rmbam = mapping.mkdup_bam
        File? uno_rmindexbam = mapping.mkdup_index

        File? mergebamfile = mergebam.mergebam
        File? mergebamindex = mergeindexstats.indexbam
        File? bkbam = merge_rmblklist.intersect_out
        File? bkindexbam = merge_bklist.indexbam
        File? rmbam = merge_markdup.mkdupbam
        File? rmindexbam = merge_mkdup.indexbam

        #MACS
        File? peakbedfile = macs.peakbedfile
        File? peakxlsfile = macs.peakxlsfile
        File? summitsfile = macs.summitsfile
        File? wigfile = macs.wigfile
        File? ctrlwigfile = macs.ctrlwigfile
        File? all_peakbedfile = all.peakbedfile
        File? all_peakxlsfile = all.peakxlsfile
        File? all_summitsfile = all.summitsfile
        File? all_wigfile = all.wigfile
        File? all_ctrlwigfile = all.ctrlwigfile
        File? nm_peakbedfile = nomodel.peakbedfile
        File? nm_peakxlsfile = nomodel.peakxlsfile
        File? nm_summitsfile = nomodel.summitsfile
        File? nm_wigfile = nomodel.wigfile
        File? nm_ctrlwigfile = nomodel.ctrlwigfile

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

        File? s_bigwig = vizsicer.bigwig
        File? s_norm_wig = vizsicer.norm_wig
        File? s_tdffile = vizsicer.tdffile

        #QC-STATS
        Array[File?]? s_qc_statsfile = s_summarystats.statsfile
        Array[File?]? s_qc_htmlfile = s_summarystats.htmlfile
        Array[File?]? s_qc_textfile = s_summarystats.textfile
        File? s_qc_mergehtml = mergehtml.mergefile

        File? s_uno_statsfile = uno_s_summarystats.statsfile
        File? s_uno_htmlfile = uno_s_summarystats.htmlfile
        File? s_uno_textfile = uno_s_summarystats.textfile

        File? uno_qc_statsfile = uno_summarystats.statsfile
        File? uno_qc_htmlfile = uno_summarystats.htmlfile
        File? uno_qc_textfile = uno_summarystats.textfile
        File? mergeqc_statsfile = merge_summarystats.statsfile
        File? mergeqc_htmlfile = merge_summarystats.htmlfile
        File? mergeqc_textfile = merge_summarystats.textfile
        File? all_summaryhtml = overallsummary.summaryhtml
        File? all_summarytxt = overallsummary.summarytxt
    }
}


