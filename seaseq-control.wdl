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
import "workflows/workflows/paired_sample_analysis.wdl" as paired

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

    # Process SRRs
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
        } # end scatter each sra

        Array[File] sample_srafile_ = flatten(fastqdump.fastqfile)
    } # end if sample_sraid

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
                    cloud=false
            }
        } # end scatter each sra

        Array[File] control_srafile_ = flatten(c_fastqdump.fastqfile)
    } # end if sample_sraid

    # Generating INDEX files
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

    # Process FASTQs
    if ( defined(sample_fastq) ) {

        Array[String] string_fastq = [1] #buffer to allow for fastq optionality
        Array[File] s_fastq = select_first([sample_fastq, string_fastq])

        Array[File] sample_fastqfile_ = s_fastq
    }

    if ( defined(control_fastq) ) {

        Array[String] c_string_fastq = [1] #buffer to allow for fastq optionality
        Array[File] c_fastq = select_first([control_fastq, c_string_fastq])

        Array[File] control_fastqfile_ = c_fastq
    }

    Array[File] bowtie_index_ = select_first([bowtie_idx_2.bowtie_indexes, bowtie_idx.bowtie_indexes, bowtie_index])
    Array[File] s_fastqfiles = flatten(select_all([sample_fastqfile_, sample_srafile_]))
    Array[File] c_fastqfiles = flatten(select_all([control_fastqfile_, control_srafile_]))

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
                    default_location='SAMPLE/' + sub(basename(eachfastq),'\.fastq\.gz|\.fq\.gz','') + '/QC/FastQC'
            }

            call util.basicfastqstats as indv_s_bfs {
                input :
                    fastqfile=eachfastq,
                    default_location='SAMPLE/' + sub(basename(eachfastq),'\.fastq\.gz|\.fq\.gz','') + '/QC/SummaryStats'
            }

            call mapping.mapping as indv_s_mapping {
                input :
                    fastqfile=eachfastq,
                    index_files=bowtie_index_,
                    metricsfile=indv_s_bfs.metrics_out,
                    blacklist=blacklist,
                    default_location='SAMPLE/' + sub(basename(eachfastq),'\.fastq\.gz|\.fq\.gz','') + '/BAM_files'
            }

            call fastqc.fastqc as indv_s_bamfqc {
                input :
                    inputfile=indv_s_mapping.sorted_bam,
                    default_location='SAMPLE/' + sub(basename(eachfastq),'\.fastq\.gz|\.fq\.gz','') + '/QC/FastQC'
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
                    default_location='SAMPLE/' + sub(basename(eachfastq),'\.fastq\.gz|\.fq\.gz','') + '/QC/SummaryStats'
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

        call samtools.index as s_mergebam_afterbklist_index_ { input: bam=s_mergebam_afterbklist }
        File s_mergebam_afterbklist_index = s_mergebam_afterbklist_index_.bam_index

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
                    default_location='CONTROL/' + sub(basename(eachfastq),'\.fastq\.gz|\.fq\.gz','') + '/QC/FastQC'
            }

            call util.basicfastqstats as indv_c_bfs {
                input :
                    fastqfile=eachfastq,
                    default_location='CONTROL/' + sub(basename(eachfastq),'\.fastq\.gz|\.fq\.gz','') + '/QC/SummaryStats'
            }

            call mapping.mapping as indv_c_mapping {
                input :
                    fastqfile=eachfastq,
                    index_files=bowtie_index_,
                    metricsfile=indv_c_bfs.metrics_out,
                    blacklist=blacklist,
                    default_location='CONTROL/' + sub(basename(eachfastq),'\.fastq\.gz|\.fq\.gz','') + '/BAM_files'
            }

            call fastqc.fastqc as indv_c_bamfqc {
                input :
                    inputfile=indv_c_mapping.sorted_bam,
                    default_location='CONTROL/' + sub(basename(eachfastq),'\.fastq\.gz|\.fq\.gz','') + '/QC/FastQC'
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
                    default_location='CONTROL/' + sub(basename(eachfastq),'\.fastq\.gz|\.fq\.gz','') + '/QC/SummaryStats'
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

        call samtools.index as c_mergebam_afterbklist_index_ { input: bam=c_mergebam_afterbklist }
        File c_mergebam_afterbklist_index = c_mergebam_afterbklist_index_.bam_index

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
                default_location='SAMPLE/' + sub(basename(s_fastqfiles[0]),'\.fastq\.gz|\.fq\.gz','') + '/QC/FastQC'
        }

        call util.basicfastqstats as uno_s_bfs {
            input :
                fastqfile=s_fastqfiles[0],
                default_location='SAMPLE/' + sub(basename(s_fastqfiles[0]),'\.fastq\.gz|\.fq\.gz','') + '/QC/SummaryStats'
        }

        call mapping.mapping as s_mapping {
            input :
                fastqfile=s_fastqfiles[0],
                index_files=bowtie_index_,
                metricsfile=uno_s_bfs.metrics_out,
                blacklist=blacklist,
                default_location='SAMPLE/' + sub(basename(s_fastqfiles[0]),'\.fastq\.gz|\.fq\.gz','') + '/BAM_files'
        }

        call fastqc.fastqc as uno_s_bamfqc {
            input :
                inputfile=s_mapping.sorted_bam,
                default_location='SAMPLE/' + sub(basename(s_fastqfiles[0]),'\.fastq\.gz|\.fq\.gz','') + '/QC/FastQC'
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
                default_location='SAMPLE/' + sub(basename(s_fastqfiles[0]),'\.fastq\.gz|\.fq\.gz','') + '/QC/SummaryStats'
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
                default_location='CONTROL/' + sub(basename(c_fastqfiles[0]),'\.fastq\.gz|\.fq\.gz','') + '/QC/FastQC'
        }

        call util.basicfastqstats as uno_c_bfs {
            input :
                fastqfile=c_fastqfiles[0],
                default_location='CONTROL/' + sub(basename(c_fastqfiles[0]),'\.fastq\.gz|\.fq\.gz','') + '/QC/SummaryStats'
        }

        call mapping.mapping as c_mapping {
            input :
                fastqfile=c_fastqfiles[0],
                index_files=bowtie_index_,
                metricsfile=uno_c_bfs.metrics_out,
                blacklist=blacklist,
                default_location='CONTROL/' + sub(basename(c_fastqfiles[0]),'\.fastq\.gz|\.fq\.gz','') + '/BAM_files'
        }

        call fastqc.fastqc as uno_c_bamfqc {
            input :
                inputfile=c_mapping.sorted_bam,
                default_location='CONTROL/' + sub(basename(c_fastqfiles[0]),'\.fastq\.gz|\.fq\.gz','') + '/QC/FastQC'
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
                default_location='CONTROL/' + sub(basename(c_fastqfiles[0]),'\.fastq\.gz|\.fq\.gz','') + '/QC/SummaryStats'
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
    # s_mergebam_afterbklist has no BAM index
    # s_mapping.bklist_bam has s_mapping.bklist_index
    # s_mapping.sorted_bam has s_mapping.bam_index
    File sample_bam = select_first([s_mergebam_afterbklist, s_mapping.bklist_bam, s_mapping.sorted_bam])
    File sample_bai = select_first([s_mergebam_afterbklist_index, s_mapping.bklist_index, s_mapping.bam_index])
    File control_bam = select_first([c_mergebam_afterbklist, c_mapping.bklist_bam, c_mapping.sorted_bam])
    File control_bai = select_first([c_mergebam_afterbklist_index, c_mapping.bklist_index, c_mapping.bam_index])

    call paired.paired_sample_analysis { input: reference=reference, blacklist=blacklist, gtf=gtf, motif_databases=motif_databases, chromsizes=samtools_faidx.chromsizes, faidx=samtools_faidx.faidx_file, sample_bam=sample_bam, sample_bai=sample_bai, control_bam=control_bam, control_bai=control_bai, results_name=results_name, run_motifs=run_motifs}

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
                bambed=paired_sample_analysis.only_s_finalbedfile,
                sppfile=paired_sample_analysis.only_s_runspp_file,
                fastqczip=select_first([uno_s_bamfqc.zipfile, string_qual]),
                bamflag=s_mapping.bam_stats,
                rmdupflag=s_mapping.mkdup_stats,
                bkflag=s_mapping.bklist_stats,
                fastqmetrics=uno_s_bfs.metrics_out,
                countsfile=paired_sample_analysis.only_s_intersectfile,
                peaksxls=paired_sample_analysis.only_s_peakxlsfile,
                outputfile = sub(basename(sample_bam),'\.sorted\.b.*$','-stats.csv'),
                outputhtml = sub(basename(sample_bam),'\.sorted\.b.*$', '-stats.html'),
                outputtext = sub(basename(sample_bam),'\.sorted\.b.*$', '-stats.txt'),
                configml = sub(basename(sample_bam),'\.sorted\.b.*$', '-config.ml')
        }

        call util.evalstats as all_summarystats {
            # SUMMARY STATISTICS of sample file (only 1 sample file provided)
            input:
                fastq_type="Comprehensive",
                bambed=paired_sample_analysis.only_s_finalbedfile,
                sppfile=paired_sample_analysis.runspp_file,
                fastqczip=select_first([uno_s_bamfqc.zipfile, string_qual]),
                bamflag=s_mapping.bam_stats,
                rmdupflag=s_mapping.mkdup_stats,
                bkflag=s_mapping.bklist_stats,
                fastqmetrics=uno_s_bfs.metrics_out,
                countsfile=paired_sample_analysis.intersect_outfile,
                peaksxls=paired_sample_analysis.peakxlsfile,
                enhancers=paired_sample_analysis.enhancers,
                superenhancers=paired_sample_analysis.super_enhancers,
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
                bambed=paired_sample_analysis.only_c_finalbedfile,
                sppfile=paired_sample_analysis.only_c_runspp_file,
                fastqczip=select_first([uno_c_bamfqc.zipfile, string_qual]),
                bamflag=c_mapping.bam_stats,
                rmdupflag=c_mapping.mkdup_stats,
                bkflag=c_mapping.bklist_stats,
                fastqmetrics=uno_c_bfs.metrics_out,
                countsfile=paired_sample_analysis.only_c_intersectfile,
                peaksxls=paired_sample_analysis.only_c_peakxlsfile,
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
                bambed=paired_sample_analysis.only_s_finalbedfile,
                sppfile=paired_sample_analysis.only_s_runspp_file,
                fastqczip=select_first([s_mergebamfqc.zipfile, string_qual]),
                bamflag=s_mergeindexstats.flagstats,
                rmdupflag=s_merge_mkdup.flagstats,
                bkflag=s_merge_bklist.flagstats,
                countsfile=paired_sample_analysis.only_s_intersectfile,
                peaksxls=paired_sample_analysis.only_s_peakxlsfile,
                outputfile = sub(basename(sample_bam),'\.sorted\.b.*$','-stats.csv'),
                outputhtml = sub(basename(sample_bam),'\.sorted\.b.*$','-stats.html'),
                outputtext = sub(basename(sample_bam),'\.sorted\.b.*$','-stats.txt'),
                configml = sub(basename(sample_bam),'\.sorted\.b.*$','-config.ml')
        }

        call util.evalstats as merge_summarystats {
            # SUMMARY STATISTICS of all samples files (more than 1 sample file provided)
            input:
                fastq_type="Comprehensive",
                bambed=paired_sample_analysis.only_s_finalbedfile,
                sppfile=paired_sample_analysis.runspp_file,
                fastqczip=select_first([s_mergebamfqc.zipfile, string_qual]),
                bamflag=s_mergeindexstats.flagstats,
                rmdupflag=s_merge_mkdup.flagstats,
                bkflag=s_merge_bklist.flagstats,
                countsfile=paired_sample_analysis.intersect_outfile,
                peaksxls=paired_sample_analysis.peakxlsfile,
                enhancers=paired_sample_analysis.enhancers,
                superenhancers=paired_sample_analysis.super_enhancers,
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
                bambed=paired_sample_analysis.only_c_finalbedfile,
                sppfile=paired_sample_analysis.only_c_runspp_file,
                fastqczip=select_first([c_mergebamfqc.zipfile, string_qual]),
                bamflag=c_mergeindexstats.flagstats,
                rmdupflag=c_merge_mkdup.flagstats,
                bkflag=c_merge_bklist.flagstats,
                countsfile=paired_sample_analysis.only_c_intersectfile,
                peaksxls=paired_sample_analysis.only_c_peakxlsfile,
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
        File? peakbedfile = paired_sample_analysis.peakbedfile
        File? peakxlsfile = paired_sample_analysis.peakxlsfile
        File? summitsfile = paired_sample_analysis.summitsfile
        File? negativexlsfile = paired_sample_analysis.negativexlsfile
        File? wigfile = paired_sample_analysis.wigfile
        File? ctrlwigfile = paired_sample_analysis.ctrlwigfile
        File? all_peakbedfile = paired_sample_analysis.all_peakbedfile
        File? all_peakxlsfile = paired_sample_analysis.all_peakxlsfile
        File? all_summitsfile = paired_sample_analysis.all_summitsfile
        File? all_negativexlsfile = paired_sample_analysis.all_negativexlsfile
        File? all_wigfile = paired_sample_analysis.all_wigfile
        File? all_ctrlwigfile = paired_sample_analysis.all_ctrlwigfile
        File? nm_peakbedfile = paired_sample_analysis.nm_peakbedfile
        File? nm_peakxlsfile = paired_sample_analysis.nm_peakxlsfile
        File? nm_summitsfile = paired_sample_analysis.nm_summitsfile
        File? nm_negativexlsfile = paired_sample_analysis.nm_negativexlsfile
        File? nm_wigfile = paired_sample_analysis.nm_wigfile
        File? nm_ctrlwigfile = paired_sample_analysis.nm_ctrlwigfile
        File? readme_peaks = paired_sample_analysis.readme_peaks

        File? only_c_peakbedfile = paired_sample_analysis.only_c_peakbedfile
        File? only_c_peakxlsfile = paired_sample_analysis.only_c_peakxlsfile
        File? only_c_summitsfile = paired_sample_analysis.only_c_summitsfile
        File? only_c_wigfile = paired_sample_analysis.only_c_wigfile
        File? only_s_peakbedfile = paired_sample_analysis.only_s_peakbedfile
        File? only_s_peakxlsfile = paired_sample_analysis.only_s_peakxlsfile
        File? only_s_summitsfile = paired_sample_analysis.only_s_summitsfile
        File? only_s_wigfile = paired_sample_analysis.only_s_wigfile

        #SICER
        File? scoreisland = paired_sample_analysis.scoreisland
        File? sicer_wigfile = paired_sample_analysis.sicer_wigfile
        File? sicer_summary = paired_sample_analysis.sicer_summary
        File? sicer_fdrisland = paired_sample_analysis.sicer_fdrisland

        #ROSE
        File? pngfile = paired_sample_analysis.pngfile
        File? mapped_union = paired_sample_analysis.mapped_union
        File? mapped_stitch = paired_sample_analysis.mapped_stitch
        File? enhancers = paired_sample_analysis.enhancers
        File? super_enhancers = paired_sample_analysis.super_enhancers
        File? gff_file = paired_sample_analysis.gff_file
        File? gff_union = paired_sample_analysis.gff_union
        File? union_enhancers = paired_sample_analysis.union_enhancers
        File? stitch_enhancers = paired_sample_analysis.stitch_enhancers
        File? e_to_g_enhancers = paired_sample_analysis.e_to_g_enhancers
        File? g_to_e_enhancers = paired_sample_analysis.g_to_e_enhancers
        File? e_to_g_super_enhancers = paired_sample_analysis.e_to_g_super_enhancers
        File? g_to_e_super_enhancers = paired_sample_analysis.g_to_e_super_enhancers

        #MOTIFS
        File? flankbedfile = paired_sample_analysis.flankbedfile

        File? ame_tsv = paired_sample_analysis.ame_tsv
        File? ame_html = paired_sample_analysis.ame_html
        File? ame_seq = paired_sample_analysis.ame_seq
        File? meme = paired_sample_analysis.meme
        File? meme_summary = paired_sample_analysis.meme_summary

        File? summit_ame_tsv = paired_sample_analysis.summit_ame_tsv
        File? summit_ame_html = paired_sample_analysis.summit_ame_html
        File? summit_ame_seq = paired_sample_analysis.summit_ame_seq
        File? summit_meme = paired_sample_analysis.summit_meme
        File? summit_meme_summary = paired_sample_analysis.summit_meme_summary

        #BAM2GFF
        File? s_matrices = paired_sample_analysis.s_matrices
        File? c_matrices = paired_sample_analysis.c_matrices
        File? densityplot = paired_sample_analysis.densityplot
        File? pdf_gene = paired_sample_analysis.pdf_gene
        File? pdf_h_gene = paired_sample_analysis.pdf_h_gene
        File? png_h_gene = paired_sample_analysis.png_h_gene
        File? pdf_promoters = paired_sample_analysis.pdf_promoters
        File? pdf_h_promoters = paired_sample_analysis.pdf_h_promoters
        File? png_h_promoters = paired_sample_analysis.png_h_promoters

        #PEAKS-ANNOTATION
        File? peak_promoters = paired_sample_analysis.peak_promoters
        File? peak_genebody = paired_sample_analysis.peak_genebody
        File? peak_window = paired_sample_analysis.peak_window
        File? peak_closest = paired_sample_analysis.peak_closest
        File? peak_comparison = paired_sample_analysis.peak_comparison
        File? gene_comparison = paired_sample_analysis.gene_comparison
        File? pdf_comparison = paired_sample_analysis.pdf_comparison

        File? all_peak_promoters = paired_sample_analysis.all_peak_promoters
        File? all_peak_genebody = paired_sample_analysis.all_peak_genebody
        File? all_peak_window = paired_sample_analysis.all_peak_window
        File? all_peak_closest = paired_sample_analysis.all_peak_closest
        File? all_peak_comparison = paired_sample_analysis.all_peak_comparison
        File? all_gene_comparison = paired_sample_analysis.all_gene_comparison
        File? all_pdf_comparison = paired_sample_analysis.all_pdf_comparison

        File? nomodel_peak_promoters = paired_sample_analysis.nomodel_peak_promoters
        File? nomodel_peak_genebody = paired_sample_analysis.nomodel_peak_genebody
        File? nomodel_peak_window = paired_sample_analysis.nomodel_peak_window
        File? nomodel_peak_closest = paired_sample_analysis.nomodel_peak_closest
        File? nomodel_peak_comparison = paired_sample_analysis.nomodel_peak_comparison
        File? nomodel_gene_comparison = paired_sample_analysis.nomodel_gene_comparison
        File? nomodel_pdf_comparison = paired_sample_analysis.nomodel_pdf_comparison

        File? sicer_peak_promoters = paired_sample_analysis.sicer_peak_promoters
        File? sicer_peak_genebody = paired_sample_analysis.sicer_peak_genebody
        File? sicer_peak_window = paired_sample_analysis.sicer_peak_window
        File? sicer_peak_closest = paired_sample_analysis.sicer_peak_closest
        File? sicer_peak_comparison = paired_sample_analysis.sicer_peak_comparison
        File? sicer_gene_comparison = paired_sample_analysis.sicer_gene_comparison
        File? sicer_pdf_comparison = paired_sample_analysis.sicer_pdf_comparison

        #VISUALIZATION
        File? bigwig = paired_sample_analysis.bigwig
        File? norm_wig = paired_sample_analysis.norm_wig
        File? tdffile = paired_sample_analysis.tdffile
        File? n_bigwig = paired_sample_analysis.n_bigwig
        File? n_norm_wig = paired_sample_analysis.n_norm_wig
        File? n_tdffile = paired_sample_analysis.n_tdffile
        File? a_bigwig = paired_sample_analysis.a_bigwig
        File? a_norm_wig = paired_sample_analysis.a_norm_wig
        File? a_tdffile = paired_sample_analysis.a_tdffile

        File? c_bigwig = paired_sample_analysis.c_bigwig
        File? c_norm_wig = paired_sample_analysis.c_norm_wig
        File? c_tdffile = paired_sample_analysis.c_tdffile
        File? c_n_bigwig = paired_sample_analysis.c_n_bigwig
        File? c_n_norm_wig = paired_sample_analysis.c_n_norm_wig
        File? c_n_tdffile = paired_sample_analysis.c_n_tdffile
        File? c_a_bigwig = paired_sample_analysis.c_a_bigwig
        File? c_a_norm_wig = paired_sample_analysis.c_a_norm_wig
        File? c_a_tdffile = paired_sample_analysis.c_a_tdffile

        File? s_bigwig = paired_sample_analysis.s_bigwig
        File? s_norm_wig = paired_sample_analysis.s_norm_wig
        File? s_tdffile = paired_sample_analysis.s_tdffile

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


