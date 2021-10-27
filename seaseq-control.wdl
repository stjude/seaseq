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
        Array[File]? bowtie_index
        Array[File]? motif_databases

        # group: input_genomic_data
        Array[String]? sample_sraid
        Array[File]? sample_fastq
        Array[String]? control_sraid
        Array[File]? control_fastq

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
            patterns: ["*.fq.gz", "*.fastq.gz"]
        }
        control_sraid: {
            description: 'One or more input/control SRA (Sequence Read Archive) run identifiers',
            group: 'input_genomic_data',
            help: 'Input publicly available FASTQs (SRRs). Multiple SRRs are separated by commas (,).',
            example: 'SRR12345678'
        }
        control_fastq: {
            description: 'One or more input/control FASTQs',
            group: 'input_genomic_data',
            help: 'Upload zipped FASTQ files.',
            patterns: ["*fq", "*.fq.gz", "*.fastq", "*.fastq.gz"]
        }
        results_name: {
            description: 'Experiment results custom name',
            group: 'analysis_parameter',
            help: 'Input preferred analysis results name.',
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
                    cloud=false
            }
        }

        Array[File] sample_fastqfile_ = flatten(fastqdump.fastqfile)
    }

    if ( defined(control_sraid) ) {
        # Download control file(s) from SRA database
        # outputs:
        #    fastqdump.fastqfile : downloaded sample files in fastq.gz format
        Array[String] c_sra = [1] #buffer to allow for sra_id optionality
        Array[String] c_sraid = select_first([control_sraid, c_sra])
        scatter (eachsra in c_sraid) {
            call sra.fastqdump as c_fastqdump{
                input :
                    sra_id=eachsra,
                    cloud="false"
            }
        }

        Array[File] control_fastqfile_ = flatten(c_fastqdump.fastqfile)
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
    Array[File] c_fastqfiles = flatten(select_all([control_fastq, control_fastqfile_]))

### ------------------------------------------------- ###
### ---------------- S E C T I O N 2 ---------------- ###
### ---- A: analysis if multiple FASTQs provided ---- ###
### ------------------------------------------------- ###

    # if multiple fastqfiles are provided
    Boolean multi_sample_fastq = if length(s_fastqfiles) > 1 then true else false
    Boolean one_sample_fastq = if length(s_fastqfiles) == 1 then true else false
    Boolean multi_control_fastq = if length(c_fastqfiles) > 1 then true else false
    Boolean one_control_fastq = if length(c_fastqfiles) == 1 then true else false

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
            
            call fastqc.fastqc as indv_s_fastqc {
                input :
                    inputfile=eachfastq,
                    default_location='SAMPLE/' + sub(basename(eachfastq),'\.f.*q\.gz','') + '/QC/FastQC'
            }

            call util.basicfastqstats as indv_s_bfs {
                input :
                    fastqfile=eachfastq,
                    default_location='SAMPLE/' + sub(basename(eachfastq),'\.f.*q\.gz','') + '/QC/SummaryStats'
            }

            call mapping.mapping as indv_s_mapping {
                input :
                    fastqfile=eachfastq,
                    index_files=bowtie_index_,
                    metricsfile=indv_s_bfs.metrics_out,
                    blacklist=blacklist,
                    default_location='SAMPLE/' + sub(basename(eachfastq),'\.f.*q\.gz','') + '/BAM_files'
            }

            call fastqc.fastqc as indv_s_bamfqc {
                input :
                    inputfile=indv_s_mapping.sorted_bam,
                    default_location='SAMPLE/' + sub(basename(eachfastq),'\.f.*q\.gz','') + '/QC/FastQC'
            }

            call runspp.runspp as indv_s_runspp {
                input:
                    bamfile=select_first([indv_s_mapping.bklist_bam, indv_s_mapping.sorted_bam])
            }

            call bedtools.bamtobed as indv_s_bamtobed {
                input:
                    bamfile=select_first([indv_s_mapping.bklist_bam, indv_s_mapping.sorted_bam])
            }

            call util.evalstats as indv_s_summarystats {
                input:
                    fastq_type="Sample FASTQ",
                    bambed=indv_s_bamtobed.bedfile,
                    sppfile=indv_s_runspp.spp_out,
                    fastqczip=indv_s_fastqc.zipfile,
                    bamflag=indv_s_mapping.bam_stats,
                    rmdupflag=indv_s_mapping.mkdup_stats,
                    bkflag=indv_s_mapping.bklist_stats,
                    fastqmetrics=indv_s_bfs.metrics_out,
                    default_location='SAMPLE/' + sub(basename(eachfastq),'\.f.*q\.gz','') + '/QC/SummaryStats'
            }
        } # end scatter (for each sample fastq)

        # MERGE BAM FILES
        # Execute analysis on merge bam file
        # Analysis executed:
        #   Merge BAM (if more than 1 fastq is provided)
        #   FastQC on Merge BAM (AllCasesMerge_<number>_mapped)

        # merge bam files and perform fasTQC if more than one is provided
        call util.mergehtml as s_mergehtml {
            input:
                htmlfiles=indv_s_summarystats.xhtml,
                txtfiles=indv_s_summarystats.textfile,
                default_location='SAMPLE',
                outputfile = 'AllCases_' + length(s_fastqfiles) + '_seaseq-summary-stats.html'
        }

        call samtools.mergebam as s_mergebam {
            input:
                bamfiles=indv_s_mapping.sorted_bam,
                default_location = 'SAMPLE/AllCasesMerge_' + length(indv_s_mapping.sorted_bam) + '_mapped' + '/BAM_files',
                outputfile = 'AllCasesMerge_' + length(s_fastqfiles) + '_mapped.sorted.bam'
        }

        call fastqc.fastqc as s_mergebamfqc {
            input:
	        inputfile=s_mergebam.mergebam,
                default_location='SAMPLE/' + sub(basename(s_mergebam.mergebam),'\.sorted\.b.*$','') + '/QC/FastQC'
        }

        call samtools.indexstats as s_mergeindexstats {
            input:
                bamfile=s_mergebam.mergebam,
                default_location='SAMPLE/' + sub(basename(s_mergebam.mergebam),'\.sorted\.b.*$','') + '/BAM_files'
        }

        if ( defined(blacklist) ) {
            # remove blacklist regions
            String string_blacklist = "" #buffer to allow for blacklist optionality
            File blacklist_ = select_first([blacklist, string_blacklist])
            call bedtools.intersect as s_merge_rmblklist {
                input :
                    fileA=s_mergebam.mergebam,
                    fileB=blacklist_,
                    default_location='SAMPLE/' + sub(basename(s_mergebam.mergebam),'\.sorted\.b.*$','') + '/BAM_files',
                    nooverlap=true
            }
            call samtools.indexstats as s_merge_bklist {
                input :
                    bamfile=s_merge_rmblklist.intersect_out,
                    default_location='SAMPLE/' + sub(basename(s_mergebam.mergebam),'\.sorted\.b.*$','') + '/BAM_files'
            }
        } # end if blacklist provided

        File s_mergebam_afterbklist = select_first([s_merge_rmblklist.intersect_out, s_mergebam.mergebam])

        call samtools.markdup as s_merge_markdup {
            input :
                bamfile=s_mergebam_afterbklist,
                default_location='SAMPLE/' + sub(basename(s_mergebam_afterbklist),'\.sorted\.b.*$','') + '/BAM_files'
        }

        call samtools.indexstats as s_merge_mkdup {
            input :
                bamfile=s_merge_markdup.mkdupbam,
                default_location='SAMPLE/' + sub(basename(s_mergebam_afterbklist),'\.sorted\.b.*$','') + '/BAM_files'
        }
    } # end if length(fastqfiles) > 1: multi_sample_fastq

    # CONTROL FASTQ files
    if ( multi_control_fastq ) {
        scatter (eachfastq in c_fastqfiles) {
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

            call fastqc.fastqc as indv_c_fastqc {
                input :
                    inputfile=eachfastq,
                    default_location='CONTROL/' + sub(basename(eachfastq),'\.f.*q\.gz','') + '/QC/FastQC'
            }

            call util.basicfastqstats as indv_c_bfs {
                input :
                    fastqfile=eachfastq,
                    default_location='CONTROL/' + sub(basename(eachfastq),'\.f.*q\.gz','') + '/QC/SummaryStats'
            }

            call mapping.mapping as indv_c_mapping {
                input :
                    fastqfile=eachfastq,
                    index_files=bowtie_index_,
                    metricsfile=indv_c_bfs.metrics_out,
                    blacklist=blacklist,
                    default_location='CONTROL/' + sub(basename(eachfastq),'\.f.*q\.gz','') + '/BAM_files'
            }

            call fastqc.fastqc as indv_c_bamfqc {
                input :
                    inputfile=indv_c_mapping.sorted_bam,
                    default_location='CONTROL/' + sub(basename(eachfastq),'\.f.*q\.gz','') + '/QC/FastQC'
            }

            call runspp.runspp as indv_c_runspp {
                input:
                    bamfile=select_first([indv_c_mapping.bklist_bam, indv_c_mapping.sorted_bam])
            }

            call bedtools.bamtobed as indv_c_bamtobed {
                input:
                    bamfile=select_first([indv_c_mapping.bklist_bam, indv_c_mapping.sorted_bam])
            }

            call util.evalstats as indv_c_summarystats {
                input:
                    fastq_type="Control FASTQ",
                    bambed=indv_c_bamtobed.bedfile,
                    sppfile=indv_c_runspp.spp_out,
                    fastqczip=indv_c_fastqc.zipfile,
                    bamflag=indv_c_mapping.bam_stats,
                    rmdupflag=indv_c_mapping.mkdup_stats,
                    bkflag=indv_c_mapping.bklist_stats,
                    fastqmetrics=indv_c_bfs.metrics_out,
                    default_location='CONTROL/' + sub(basename(eachfastq),'\.f.*q\.gz','') + '/QC/SummaryStats'
            }
        } # end scatter (for each control fastq)

        # MERGE BAM FILES
        # Execute analysis on merge bam file
        # Analysis executed:
        #   Merge BAM (if more than 1 fastq is provided)
        #   FastQC on Merge BAM (AllCtrlsMerge_<number>_mapped)

        # merge bam files and perform fasTQC if more than one is provided
        call util.mergehtml as c_mergehtml {
            input:
                fastq_type="Control FASTQs",
                htmlfiles=indv_c_summarystats.xhtml,
                txtfiles=indv_c_summarystats.textfile,
                default_location='CONTROL',
                outputfile = 'AllCtrls_' + length(c_fastqfiles) + '_seaseq-summary-stats.html'
        }

        call samtools.mergebam as c_mergebam {
            input:
                bamfiles=indv_c_mapping.sorted_bam,
                default_location = 'CONTROL/' + 'AllCtrlsMerge_' + length(indv_c_mapping.sorted_bam) + '_mapped' + '/BAM_files',
                outputfile = 'AllCtrlsMerge_' + length(c_fastqfiles) + '_mapped.sorted.bam'
        }

        call fastqc.fastqc as c_mergebamfqc {
            input:
	        inputfile=c_mergebam.mergebam,
                default_location='CONTROL/' + sub(basename(c_mergebam.mergebam),'\.sorted\.b.*$','') + '/QC/FastQC'
        }

        call samtools.indexstats as c_mergeindexstats {
            input:
                bamfile=c_mergebam.mergebam,
                default_location='CONTROL/' + sub(basename(c_mergebam.mergebam),'\.sorted\.b.*$','') + '/BAM_files'
        }

        if ( defined(blacklist) ) {
            # remove blacklist regions
            String c_string_blacklist = "" #buffer to allow for blacklist optionality
            File c_blacklist_ = select_first([blacklist, c_string_blacklist])
            call bedtools.intersect as c_merge_rmblklist {
                input :
                    fileA=c_mergebam.mergebam,
                    fileB=c_blacklist_,
                    default_location='CONTROL/' + sub(basename(c_mergebam.mergebam),'\.sorted\.b.*$','') + '/BAM_files',
                    nooverlap=true
            }
            call samtools.indexstats as c_merge_bklist {
                input :
                    bamfile=c_merge_rmblklist.intersect_out,
                    default_location='CONTROL/' + sub(basename(c_mergebam.mergebam),'\.sorted\.b.*$','') + '/BAM_files'
            }
        } # end if blacklist provided

        File c_mergebam_afterbklist = select_first([c_merge_rmblklist.intersect_out, c_mergebam.mergebam])

        call samtools.markdup as c_merge_markdup {
            input :
                bamfile=c_mergebam_afterbklist,
                default_location='CONTROL/' + sub(basename(c_mergebam_afterbklist),'\.sorted\.b.*$','') + '/BAM_files'
        }

        call samtools.indexstats as c_merge_mkdup {
            input :
                bamfile=c_merge_markdup.mkdupbam,
                default_location='CONTROL/' + sub(basename(c_mergebam_afterbklist),'\.sorted\.b.*$','') + '/BAM_files'
        }
    } # end if length(fastqfiles) > 1: multi_control_fastq


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

        call fastqc.fastqc as uno_s_fastqc {
            input :
                inputfile=s_fastqfiles[0],
                default_location='SAMPLE/' + sub(basename(s_fastqfiles[0]),'\.f.*q\.gz','') + '/QC/FastQC'
        }

        call util.basicfastqstats as uno_s_bfs {
            input :
                fastqfile=s_fastqfiles[0],
                default_location='SAMPLE/' + sub(basename(s_fastqfiles[0]),'\.f.*q\.gz','') + '/QC/SummaryStats'
        }

        call mapping.mapping as s_mapping {
            input :
                fastqfile=s_fastqfiles[0],
                index_files=bowtie_index_,
                metricsfile=uno_s_bfs.metrics_out,
                blacklist=blacklist,
                default_location='SAMPLE/' + sub(basename(s_fastqfiles[0]),'\.f.*q\.gz','') + '/BAM_files'
        }

        call fastqc.fastqc as uno_s_bamfqc {
            input :
                inputfile=s_mapping.sorted_bam,
                default_location='SAMPLE/' + sub(basename(s_fastqfiles[0]),'\.f.*q\.gz','') + '/QC/FastQC'
        }

        call runspp.runspp as uno_s_runspp {
            input:
                bamfile=select_first([s_mapping.bklist_bam, s_mapping.sorted_bam])
        }

        call bedtools.bamtobed as uno_s_bamtobed {
            input:
                bamfile=select_first([s_mapping.bklist_bam, s_mapping.sorted_bam])
        }

        call util.evalstats as uno_s_summarystats {
            input:
                fastq_type="Sample FASTQ",
                bambed=uno_s_bamtobed.bedfile,
                sppfile=uno_s_runspp.spp_out,
                fastqczip=uno_s_fastqc.zipfile,
                bamflag=s_mapping.bam_stats,
                rmdupflag=s_mapping.mkdup_stats,
                bkflag=s_mapping.bklist_stats,
                fastqmetrics=uno_s_bfs.metrics_out,
                default_location='SAMPLE/' + sub(basename(s_fastqfiles[0]),'\.f.*q\.gz','') + '/QC/SummaryStats'
        }
    } # end if length(fastqfiles) == 1: one_sample_fastq

    # CONTROL FASTQ file
    if ( one_control_fastq ) {
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
        
        call fastqc.fastqc as uno_c_fastqc {
            input :
                inputfile=c_fastqfiles[0],
                default_location='CONTROL/' + sub(basename(c_fastqfiles[0]),'\.f.*q\.gz','') + '/QC/FastQC'
        }

        call util.basicfastqstats as uno_c_bfs {
            input :
                fastqfile=c_fastqfiles[0],
                default_location='CONTROL/' + sub(basename(c_fastqfiles[0]),'\.f.*q\.gz','') + '/QC/SummaryStats'
        }

        call mapping.mapping as c_mapping {
            input :
                fastqfile=c_fastqfiles[0],
                index_files=bowtie_index_,
                metricsfile=uno_c_bfs.metrics_out,
                blacklist=blacklist,
                default_location='CONTROL/' + sub(basename(c_fastqfiles[0]),'\.f.*q\.gz','') + '/BAM_files'
        }

        call fastqc.fastqc as uno_c_bamfqc {
            input :
                inputfile=c_mapping.sorted_bam,
                default_location='CONTROL/' + sub(basename(c_fastqfiles[0]),'\.f.*q\.gz','') + '/QC/FastQC'
        }

        call runspp.runspp as uno_c_runspp {
            input:
                bamfile=select_first([c_mapping.bklist_bam, c_mapping.sorted_bam])
        }

        call bedtools.bamtobed as uno_c_bamtobed {
            input:
                bamfile=select_first([c_mapping.bklist_bam, c_mapping.sorted_bam])
        }

        call util.evalstats as uno_c_summarystats {
            input:
                fastq_type="Control FASTQ",
                bambed=uno_c_bamtobed.bedfile,
                sppfile=uno_c_runspp.spp_out,
                fastqczip=uno_c_fastqc.zipfile,
                bamflag=c_mapping.bam_stats,
                rmdupflag=c_mapping.mkdup_stats,
                bkflag=c_mapping.bklist_stats,
                fastqmetrics=uno_c_bfs.metrics_out,
                default_location='CONTROL/' + sub(basename(c_fastqfiles[0]),'\.f.*q\.gz','') + '/QC/SummaryStats'
        }
    } # end if length(fastqfiles) == 1: one_control_fastq

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
    File sample_bam = select_first([s_mergebam_afterbklist, s_mapping.bklist_bam, s_mapping.sorted_bam])
    File control_bam = select_first([c_mergebam_afterbklist, c_mapping.bklist_bam, c_mapping.sorted_bam])

    call macs.macs {
        input :
            bamfile=sample_bam,
            control=control_bam,
            pvalue = "1e-9",
            keep_dup="auto",
            output_name = if defined(results_name) then results_name + '-p9_kd-auto' else basename(sample_bam,'\.bam') + '+control-p9_kd-auto',
            default_location = if defined(results_name) then results_name + '/PEAKS/NARROW_peaks/' + results_name + '-p9_kd-auto' else sub(basename(sample_bam),'\.sorted\.b.*$','') + '+control/PEAKS/NARROW_peaks/' + basename(sample_bam,'\.bam') + '+control-p9_kd-auto'
    }

    call util.addreadme {
        input :
            default_location = if defined(results_name) then results_name + '/PEAKS' else sub(basename(sample_bam),'\.sorted\.b.*$','') + '+control/PEAKS/'
    }
    
    call macs.macs as all {
        input :
            bamfile=sample_bam,
            control=control_bam,
            pvalue = "1e-9",
            keep_dup="all",
            output_name = if defined(results_name) then results_name + '-p9_kd-all' else basename(sample_bam,'\.bam') + '+control-p9_kd-all',
            default_location = if defined(results_name) then results_name + '/PEAKS/NARROW_peaks/' + results_name + '-p9_kd-all' else sub(basename(sample_bam),'\.sorted\.b.*$','') + '+control/PEAKS/NARROW_peaks/' + basename(sample_bam,'\.bam') + '+control-p9_kd-all'
    }

    call macs.macs as nomodel {
        input :
            bamfile=sample_bam,
            control=control_bam,
            nomodel=true,
            output_name = if defined(results_name) then results_name + '-nm' else basename(sample_bam,'\.bam') + '+control-nm',
            default_location = if defined(results_name) then results_name + '/PEAKS/NARROW_peaks/' + results_name + '-nm' else sub(basename(sample_bam),'\.sorted\.b.*$','') + '+control/PEAKS/NARROW_peaks/' + basename(sample_bam,'\.bam') + '+control-nm'
    }

    call bamtogff.bamtogff {
        input :
            gtffile=gtf,
            chromsizes=samtools_faidx.chromsizes,
            bamfile=select_first([s_merge_markdup.mkdupbam, s_mapping.mkdup_bam]),
            bamindex=select_first([s_merge_mkdup.indexbam, s_mapping.mkdup_index]),
            control_bamfile=select_first([c_merge_markdup.mkdupbam, c_mapping.mkdup_bam]),
            control_bamindex=select_first([c_merge_mkdup.indexbam, c_mapping.mkdup_index]),
            samplename=if defined(results_name) then results_name else basename(sample_bam,'\.bam') + '+control',
            default_location=if defined(results_name) then results_name + '/BAM_Density' else sub(basename(sample_bam),'\.sorted\.b.*$','') + '+control/BAM_Density'
    }

    call bedtools.bamtobed as s_forsicerbed {
        input :
            bamfile=select_first([s_merge_markdup.mkdupbam, s_mapping.mkdup_bam])
    }

    call bedtools.bamtobed as c_forsicerbed {
        input :
            bamfile=select_first([c_merge_markdup.mkdupbam, c_mapping.mkdup_bam])
    }

    call sicer.sicer {
        input :
            bedfile=s_forsicerbed.bedfile,
            control_bed=c_forsicerbed.bedfile,
            chromsizes=samtools_faidx.chromsizes,
            outputname=if defined(results_name) then results_name else basename(s_forsicerbed.bedfile,'\.bed') + '+control',
            default_location=if defined(results_name) then results_name + '/PEAKS/BROAD_peaks' else sub(basename(sample_bam),'\.sorted\.b.*$','') + '+control/PEAKS/BROAD_peaks'
    }

    call rose.rose {
        input :
            gtffile=gtf,
            bamfile=sample_bam,
            bamindex=select_first([s_merge_bklist.indexbam, s_mergeindexstats.indexbam, s_mapping.bklist_index, s_mapping.bam_index]),
            control=control_bam,
            controlindex=select_first([c_merge_bklist.indexbam, c_mergeindexstats.indexbam, c_mapping.bklist_index, c_mapping.bam_index]),
            bedfile_auto=macs.peakbedfile,
            bedfile_all=all.peakbedfile,
            default_location=if defined(results_name) then results_name + '/PEAKS/STITCHED_peaks' else sub(basename(sample_bam),'\.sorted\.b.*$','') + '+control/PEAKS/STITCHED_peaks'
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
            chromsizes=samtools_faidx.chromsizes,
            control=true,
            xlsfile=macs.peakxlsfile,
            default_location=if defined(results_name) then results_name + '/COVERAGE_files/NARROW_peaks/' + sub(basename(macs.peakbedfile),'\_peaks.bed','') + '/control' else sub(basename(sample_bam),'\.sorted\.b.*$','') + '+control/COVERAGE_files/NARROW_peaks/' + sub(basename(macs.peakbedfile),'\_peaks.bed','') + '/control'
    }

    call viz.visualization as c_vizall {
        input:
            wigfile=select_first([all.ctrlwigfile, string_ctrlwig]),
            chromsizes=samtools_faidx.chromsizes,
            control=true,
            xlsfile=all.peakxlsfile,
            default_location=if defined(results_name) then results_name + '/COVERAGE_files/NARROW_peaks/' + sub(basename(all.peakbedfile),'\_peaks.bed','') + '/control' else sub(basename(sample_bam),'\.sorted\.b.*$','') + '+control/COVERAGE_files/NARROW_peaks/' + sub(basename(all.peakbedfile),'\_peaks.bed','') + '/control'
    }
    call viz.visualization as c_viznomodel {
        input:
            wigfile=select_first([nomodel.ctrlwigfile, string_ctrlwig]),
            chromsizes=samtools_faidx.chromsizes,
            control=true,
            xlsfile=nomodel.peakxlsfile,
            default_location=if defined(results_name) then results_name + '/COVERAGE_files/NARROW_peaks/' + sub(basename(nomodel.peakbedfile),'\_peaks.bed','') + '/control' else sub(basename(sample_bam),'\.sorted\.b.*$','') + '+control/COVERAGE_files/NARROW_peaks/' + sub(basename(nomodel.peakbedfile),'\_peaks.bed','') + '/control'
    }

    call util.peaksanno {
        input :
            gtffile=gtf,
            bedfile=macs.peakbedfile,
            chromsizes=samtools_faidx.chromsizes,
            summitfile=macs.summitsfile,
            default_location=if defined(results_name) then results_name + '/PEAKS_Annotation/NARROW_peaks/' + sub(basename(macs.peakbedfile),'\_peaks.bed','') else sub(basename(sample_bam),'\.sorted\.b.*$','') + '+control/PEAKS_Annotation/NARROW_peaks/' + sub(basename(macs.peakbedfile),'\_peaks.bed','')
    }

    call util.peaksanno as all_peaksanno {
        input :
            gtffile=gtf,
            bedfile=all.peakbedfile,
            chromsizes=samtools_faidx.chromsizes,
            summitfile=all.summitsfile,
            default_location=if defined(results_name) then results_name + '/PEAKS_Annotation/NARROW_peaks/' + sub(basename(all.peakbedfile),'\_peaks.bed','') else sub(basename(sample_bam),'\.sorted\.b.*$','') + '+control/PEAKS_Annotation/NARROW_peaks/' + sub(basename(all.peakbedfile),'\_peaks.bed','')
    }

    call util.peaksanno as nomodel_peaksanno {
        input :
            gtffile=gtf,
            bedfile=nomodel.peakbedfile,
            chromsizes=samtools_faidx.chromsizes,
            summitfile=nomodel.summitsfile,
            default_location=if defined(results_name) then results_name + '/PEAKS_Annotation/NARROW_peaks/' + sub(basename(nomodel.peakbedfile),'\_peaks.bed','') else sub(basename(sample_bam),'\.sorted\.b.*$','') + '+control/PEAKS_Annotation/NARROW_peaks/' + sub(basename(nomodel.peakbedfile),'\_peaks.bed','')
    }

    call util.peaksanno as sicer_peaksanno {
        input :
            gtffile=gtf,
            bedfile=select_first([sicer.fdrisland, string_ctrlwig]),
            chromsizes=samtools_faidx.chromsizes,
            default_location=if defined(results_name) then results_name + '/PEAKS_Annotation/BROAD_peaks' else sub(basename(sample_bam),'\.sorted\.b.*$','') + '+control/PEAKS_Annotation/BROAD_peaks'
    }

    # Motif Analysis
    if (run_motifs) {
        call motifs.motifs {
            input:
                reference=reference,
                reference_index=samtools_faidx.faidx_file,
                bedfile=macs.peakbedfile,
                motif_databases=motif_databases,
                default_location=if defined(results_name) then results_name + '/MOTIFS' else sub(basename(sample_bam),'\.sorted\.b.*$','') + '+control/MOTIFS'
        }

        call util.flankbed {
            input :
                bedfile=macs.summitsfile,
                default_location=if defined(results_name) then results_name + '/MOTIFS' else sub(basename(sample_bam),'\.sorted\.b.*$','') + '+control/MOTIFS'
        }

        call motifs.motifs as flank {
            input:
                reference=reference,
                reference_index=samtools_faidx.faidx_file,
                bedfile=flankbed.flankbedfile,
                motif_databases=motif_databases,
                default_location=if defined(results_name) then results_name + '/MOTIFS' else sub(basename(sample_bam),'\.sorted\.b.*$','') + '+control/MOTIFS'
        }
    }

    call viz.visualization {
        input:
            wigfile=macs.wigfile,
            chromsizes=samtools_faidx.chromsizes,
            xlsfile=macs.peakxlsfile,
            default_location=if defined(results_name) then results_name + '/COVERAGE_files/NARROW_peaks/' + sub(basename(macs.peakbedfile),'\_peaks.bed','') else sub(basename(sample_bam),'\.sorted\.b.*$','') + '+control/COVERAGE_files/NARROW_peaks/' + sub(basename(macs.peakbedfile),'\_peaks.bed','')
    }

    call viz.visualization as vizall {
        input:
            wigfile=all.wigfile,
            chromsizes=samtools_faidx.chromsizes,
            xlsfile=all.peakxlsfile,
            default_location=if defined(results_name) then results_name + '/COVERAGE_files/NARROW_peaks/' + sub(basename(all.peakbedfile),'\_peaks.bed','') else sub(basename(sample_bam),'\.sorted\.b.*$','') + '+control/COVERAGE_files/NARROW_peaks/' + sub(basename(all.peakbedfile),'\_peaks.bed','')
    }

    call viz.visualization as viznomodel {
        input:
            wigfile=nomodel.wigfile,
            chromsizes=samtools_faidx.chromsizes,
            xlsfile=nomodel.peakxlsfile,
            default_location=if defined(results_name) then results_name + '/COVERAGE_files/NARROW_peaks/' + sub(basename(nomodel.peakbedfile),'\_peaks.bed','') else sub(basename(sample_bam),'\.sorted\.b.*$','') + '+control/COVERAGE_files/NARROW_peaks/' + sub(basename(nomodel.peakbedfile),'\_peaks.bed','')
    }

    call viz.visualization as vizsicer {
        input:
            wigfile=sicer.wigfile,
            chromsizes=samtools_faidx.chromsizes,
            default_location=if defined(results_name) then results_name + '/COVERAGE_files/BROAD_peaks' else sub(basename(sample_bam),'\.sorted\.b.*$','') + '+control/COVERAGE_files/BROAD_peaks'
    }

    #Peak Calling for Sample BAM only
    call macs.macs as only_s_macs {
        input :
            bamfile=sample_bam,
            pvalue = "1e-9",
            keep_dup="auto",
            default_location='SAMPLE/' + sub(basename(sample_bam),'\.sorted\.b.*$','') + '/PEAKS_forQC/' + basename(sample_bam,'\.bam') + '-p9_kd-auto'
    }

    #Peak Calling for Control BAM only
    call macs.macs as only_c_macs {
        input :
            bamfile=control_bam,
            pvalue = "1e-9",
            keep_dup="auto",
            default_location='CONTROL/' + sub(basename(control_bam),'\.sorted\.b.*$','') + '/PEAKS_forQC/' + basename(control_bam,'\.bam') + '-p9_kd-auto'
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

### ---------------------------------------- ###
### ------------ S E C T I O N 4 ----------- ###
### ---------- Summary Statistics ---------- ###
### ---------------------------------------- ###

    String string_qual = "" #buffer to allow for optionality in if statement

    #SUMMARY STATISTICS
    if ( one_sample_fastq ) {
        call util.evalstats as all_s_summarystats {
            # SUMMARY STATISTICS of sample file (only 1 sample file provided)
            input:
                fastq_type="Sample",
                bambed=only_s_finalbed.bedfile,
                sppfile=only_s_runspp.spp_out,
                fastqczip=select_first([uno_s_bamfqc.zipfile, string_qual]),
                bamflag=s_mapping.bam_stats,
                rmdupflag=s_mapping.mkdup_stats,
                bkflag=s_mapping.bklist_stats,
                fastqmetrics=uno_s_bfs.metrics_out,
                countsfile=only_s_intersect.intersect_out,
                peaksxls=only_s_macs.peakxlsfile,
                outputfile = sub(basename(sample_bam),'\.sorted\.b.*$','-stats.csv'),
                outputhtml = sub(basename(sample_bam),'\.sorted\.b.*$', '-stats.html'),
                outputtext = sub(basename(sample_bam),'\.sorted\.b.*$', '-stats.txt'),
                configml = sub(basename(sample_bam),'\.sorted\.b.*$', '-config.ml')
        }

        call util.evalstats as all_summarystats {
            # SUMMARY STATISTICS of sample file (only 1 sample file provided)
            input:
                fastq_type="Comprehensive",
                bambed=only_s_finalbed.bedfile,
                sppfile=runspp.spp_out,
                fastqczip=select_first([uno_s_bamfqc.zipfile, string_qual]),
                bamflag=s_mapping.bam_stats,
                rmdupflag=s_mapping.mkdup_stats,
                bkflag=s_mapping.bklist_stats,
                fastqmetrics=uno_s_bfs.metrics_out,
                countsfile=only_s_intersect.intersect_out,
                peaksxls=macs.peakxlsfile,
                enhancers=rose.enhancers,
                superenhancers=rose.super_enhancers,
                outputfile = sub(basename(sample_bam),'\.sorted\.b.*$', '-stats.csv'),
                outputhtml = sub(basename(sample_bam),'\.sorted\.b.*$','-stats.html'),
                outputtext = sub(basename(sample_bam),'\.sorted\.b.*$', '-stats.txt'),
                configml = sub(basename(sample_bam),'\.sorted\.b.*$', '-config.ml')
        }
    } # end if one_fastq

    if ( one_control_fastq ) {
        call util.evalstats as all_c_summarystats {
            # SUMMARY STATISTICS of sample file (only 1 sample file provided)
            input:
                fastq_type="Control",
                bambed=only_c_finalbed.bedfile,
                sppfile=only_c_runspp.spp_out,
                fastqczip=select_first([uno_c_bamfqc.zipfile, string_qual]),
                bamflag=c_mapping.bam_stats,
                rmdupflag=c_mapping.mkdup_stats,
                bkflag=c_mapping.bklist_stats,
                fastqmetrics=uno_c_bfs.metrics_out,
                countsfile=only_c_intersect.intersect_out,
                peaksxls=only_c_macs.peakxlsfile,
                outputfile = sub(basename(control_bam),'\.sorted\.b.*$','-stats.csv'),
                outputhtml = sub(basename(control_bam),'\.sorted\.b.*$','-stats.html'),
                outputtext = sub(basename(control_bam),'\.sorted\.b.*$','-stats.txt'),
                configml = sub(basename(control_bam),'\.sorted\.b.*$','-config.ml')
        }
    } # end if one_fastq

    if ( multi_sample_fastq ) {
        call util.evalstats as merge_s_summarystats {
            # SUMMARY STATISTICS of all samples files (more than 1 sample file provided)
            input:
                fastq_type="Sample",
                bambed=only_s_finalbed.bedfile,
                sppfile=only_s_runspp.spp_out,
                fastqczip=select_first([s_mergebamfqc.zipfile, string_qual]),
                bamflag=s_mergeindexstats.flagstats,
                rmdupflag=s_merge_mkdup.flagstats,
                bkflag=s_merge_bklist.flagstats,
                countsfile=only_s_intersect.intersect_out,
                peaksxls=only_s_macs.peakxlsfile,
                outputfile = sub(basename(sample_bam),'\.sorted\.b.*$','-stats.csv'),
                outputhtml = sub(basename(sample_bam),'\.sorted\.b.*$','-stats.html'),
                outputtext = sub(basename(sample_bam),'\.sorted\.b.*$','-stats.txt'),
                configml = sub(basename(sample_bam),'\.sorted\.b.*$','-config.ml')
        }

        call util.evalstats as merge_summarystats {
            # SUMMARY STATISTICS of all samples files (more than 1 sample file provided)
            input:
                fastq_type="Comprehensive",
                bambed=only_s_finalbed.bedfile,
                sppfile=runspp.spp_out,
                fastqczip=select_first([s_mergebamfqc.zipfile, string_qual]),
                bamflag=s_mergeindexstats.flagstats,
                rmdupflag=s_merge_mkdup.flagstats,
                bkflag=s_merge_bklist.flagstats,
                countsfile=only_s_intersect.intersect_out,
                peaksxls=macs.peakxlsfile,
                enhancers=rose.enhancers,
                superenhancers=rose.super_enhancers,
                outputfile = sub(basename(sample_bam),'\.sorted\.b.*$','-stats.csv'),
                outputhtml = sub(basename(sample_bam),'\.sorted\.b.*$','-stats.html'),
                outputtext = sub(basename(sample_bam),'\.sorted\.b.*$','-stats.txt'),
                configml = sub(basename(sample_bam),'\.sorted\.b.*$','-config.ml')
        }
    } # end if multi_fastq

    if ( multi_control_fastq ) {
        call util.evalstats as merge_c_summarystats {
            # SUMMARY STATISTICS of all samples files (more than 1 sample file provided)
            input:
                fastq_type="Control",
                bambed=only_c_finalbed.bedfile,
                sppfile=only_c_runspp.spp_out,
                fastqczip=select_first([c_mergebamfqc.zipfile, string_qual]),
                bamflag=c_mergeindexstats.flagstats,
                rmdupflag=c_merge_mkdup.flagstats,
                bkflag=c_merge_bklist.flagstats,
                countsfile=only_c_intersect.intersect_out,
                peaksxls=only_c_macs.peakxlsfile,
                outputfile = sub(basename(control_bam),'\.sorted\.b.*$','-stats.csv'),
                outputhtml = sub(basename(control_bam),'\.sorted\.b.*$','-stats.html'),
                outputtext = sub(basename(control_bam),'\.sorted\.b.*$','-stats.txt'),
                configml = sub(basename(control_bam),'\.sorted\.b.*$','-config.ml')
        }
    } # end if multi_fastq

    call util.concatstats {
        # SUMMARY STATISTICS of all samples files (more than 1 sample file provided)
        input:
            sample_config=select_first([all_s_summarystats.configfile, merge_s_summarystats.configfile]),
            control_config=select_first([all_c_summarystats.configfile, merge_c_summarystats.configfile]),
            overall_config=select_first([all_summarystats.configfile, merge_summarystats.configfile]),
            outputfile=if defined(results_name) then results_name else sub(basename(sample_bam),'\.sorted\.b.*$','') + '+control',
            default_location=if defined(results_name) then results_name + '/QC/SummaryStats' else sub(basename(sample_bam),'\.sorted\.b.*$','') + '+control/QC/SummaryStats'
    }

    call util.summaryreport as overallsummary {
        # Presenting all quality stats for the analysis
        input:
            controlqc_html=select_first([uno_c_summarystats.xhtml, c_mergehtml.xhtml]),
            sampleqc_html=select_first([uno_s_summarystats.xhtml, s_mergehtml.xhtml]),
            overallqc_html=concatstats.xhtml,
            controlqc_txt=select_first([uno_c_summarystats.textfile, c_mergehtml.mergetxt, string_qual]),
            sampleqc_txt=select_first([uno_s_summarystats.textfile, s_mergehtml.mergetxt]),
            overallqc_txt=concatstats.textfile,
    }

    output {
        #FASTQC
        Array[File?]? indv_s_htmlfile = indv_s_fastqc.htmlfile
        Array[File?]? indv_s_zipfile = indv_s_fastqc.zipfile
        Array[File?]? indv_s_bam_htmlfile = indv_s_bamfqc.htmlfile
        Array[File?]? indv_s_bam_zipfile = indv_s_bamfqc.zipfile
        Array[File?]? indv_c_htmlfile = indv_c_fastqc.htmlfile
        Array[File?]? indv_c_zipfile = indv_c_fastqc.zipfile
        Array[File?]? indv_c_bam_htmlfile = indv_c_bamfqc.htmlfile
        Array[File?]? indv_c_bam_zipfile = indv_c_bamfqc.zipfile

        File? s_mergebam_htmlfile = s_mergebamfqc.htmlfile
        File? s_mergebam_zipfile = s_mergebamfqc.zipfile
        File? c_mergebam_htmlfile = c_mergebamfqc.htmlfile
        File? c_mergebam_zipfile = c_mergebamfqc.zipfile

        File? uno_s_htmlfile = uno_s_fastqc.htmlfile
        File? uno_s_zipfile = uno_s_fastqc.zipfile
        File? uno_s_bam_htmlfile = uno_s_bamfqc.htmlfile
        File? uno_s_bam_zipfile = uno_s_bamfqc.zipfile
        File? uno_c_htmlfile = uno_c_fastqc.htmlfile
        File? uno_c_zipfile = uno_c_fastqc.zipfile
        File? uno_c_bam_htmlfile = uno_c_bamfqc.htmlfile
        File? uno_c_bam_zipfile = uno_c_bamfqc.zipfile

        #BASICMETRICS
        Array[File?]? s_metrics_out = indv_s_bfs.metrics_out
        File? uno_s_metrics_out = uno_s_bfs.metrics_out
        Array[File?]? c_metrics_out = indv_c_bfs.metrics_out
        File? uno_c_metrics_out = uno_c_bfs.metrics_out

        #BAMFILES
        Array[File?]? indv_s_sortedbam = indv_s_mapping.sorted_bam
        Array[File?]? indv_s_indexbam = indv_s_mapping.bam_index
        Array[File?]? indv_s_bkbam = indv_s_mapping.bklist_bam
        Array[File?]? indv_s_bkindexbam = indv_s_mapping.bklist_index
        Array[File?]? indv_s_rmbam = indv_s_mapping.mkdup_bam
        Array[File?]? indv_s_rmindexbam = indv_s_mapping.mkdup_index
        Array[File?]? indv_c_sortedbam = indv_c_mapping.sorted_bam
        Array[File?]? indv_c_indexbam = indv_c_mapping.bam_index
        Array[File?]? indv_c_bkbam = indv_c_mapping.bklist_bam
        Array[File?]? indv_c_bkindexbam = indv_c_mapping.bklist_index
        Array[File?]? indv_c_rmbam = indv_c_mapping.mkdup_bam
        Array[File?]? indv_c_rmindexbam = indv_c_mapping.mkdup_index

        File? uno_s_sortedbam = s_mapping.sorted_bam
        File? uno_s_indexstatsbam = s_mapping.bam_index
        File? uno_s_bkbam = s_mapping.bklist_bam
        File? uno_s_bkindexbam = s_mapping.bklist_index
        File? uno_s_rmbam = s_mapping.mkdup_bam
        File? uno_s_rmindexbam = s_mapping.mkdup_index
        File? uno_c_sortedbam = c_mapping.sorted_bam
        File? uno_c_indexstatsbam = c_mapping.bam_index
        File? uno_c_bkbam = c_mapping.bklist_bam
        File? uno_c_bkindexbam = c_mapping.bklist_index
        File? uno_c_rmbam = c_mapping.mkdup_bam
        File? uno_c_rmindexbam = c_mapping.mkdup_index

        File? s_mergebamfile = s_mergebam.mergebam
        File? s_mergebamindex = s_mergeindexstats.indexbam
        File? s_bkbam = s_merge_rmblklist.intersect_out
        File? s_bkindexbam = s_merge_bklist.indexbam
        File? s_rmbam = s_merge_markdup.mkdupbam
        File? s_rmindexbam = s_merge_mkdup.indexbam
        File? c_mergebamfile = c_mergebam.mergebam
        File? c_mergebamindex = c_mergeindexstats.indexbam
        File? c_bkbam = c_merge_rmblklist.intersect_out
        File? c_bkindexbam = c_merge_bklist.indexbam
        File? c_rmbam = c_merge_markdup.mkdupbam
        File? c_rmindexbam = c_merge_mkdup.indexbam
        
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

        #QC-STATS
        Array[File?]? s_qc_statsfile = indv_s_summarystats.statsfile
        Array[File?]? s_qc_htmlfile = indv_s_summarystats.htmlfile
        Array[File?]? s_qc_textfile = indv_s_summarystats.textfile
        File? s_qc_mergehtml = s_mergehtml.mergefile
        Array[File?]? c_qc_statsfile = indv_c_summarystats.statsfile
        Array[File?]? c_qc_htmlfile = indv_c_summarystats.htmlfile
        Array[File?]? c_qc_textfile = indv_c_summarystats.textfile
        File? c_qc_mergehtml = c_mergehtml.mergefile

        File? s_uno_statsfile = uno_s_summarystats.statsfile
        File? s_uno_htmlfile = uno_s_summarystats.htmlfile
        File? s_uno_textfile = uno_s_summarystats.textfile
        File? c_uno_statsfile = uno_c_summarystats.statsfile
        File? c_uno_htmlfile = uno_c_summarystats.htmlfile
        File? c_uno_textfile = uno_c_summarystats.textfile

        File? statsfile = concatstats.statsfile
        File? htmlfile = concatstats.htmlfile
        File? textfile = concatstats.textfile

        File? summaryhtml = overallsummary.summaryhtml
        File? summarytxt = overallsummary.summarytxt
    }
}


