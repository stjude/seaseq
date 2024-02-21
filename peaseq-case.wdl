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
import "workflows/tasks/seaseq_util.wdl" as util
import "workflows/workflows/visualization.wdl" as viz
import "workflows/workflows/mapping.wdl"
import "workflows/tasks/runspp.wdl"
import "workflows/tasks/sortbed.wdl"
import "workflows/tasks/sratoolkit.wdl" as sra
import "workflows/tasks/peaseq_util.wdl"

workflow peaseq {
    String pipeline_ver = 'v1.0.0'

    meta {
        title: 'PEAseq Analysis'
        summary: 'Paired-End Antibody Sequencing (PEAseq) Pipeline'
        description: 'A comprehensive automated computational pipeline for all ChIP-Seq/CUT&RUN data analysis.'
        version: '1.0.0'
        details: {
            citation: 'https://doi.org/10.1186/s12859-022-04588-z',
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
        Array[File]? sample_R1_fastq
        Array[File]? sample_R2_fastq

        # group: analysis_parameter
        # Read Mapping parameters for bowtie
        Int? insertsize = 600
        String? strandedness = "fr"

        # Additional options
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
        sample_R1_fastq: {
            description: 'One or more sample R1 FASTQs',
            group: 'input_genomic_data',
            help: 'Upload zipped FASTQ files.',
            patterns: ["*.fq.gz", "*.fastq.gz"]
        }
        sample_R2_fastq: {
            description: 'One or more sample R2 FASTQs',
            group: 'input_genomic_data',
            help: 'Upload zipped FASTQ files.',
            patterns: ["*.fq.gz", "*.fastq.gz"]
        }
        results_name: {
            description: 'Experiment results custom name',
            group: 'analysis_parameter',
            help: 'Input preferred analysis results name (recommended if multiple FASTQs are provided).',
            example: 'AllMerge_mapped'
        }
        run_motifs: {
            description: 'Perform Motif Analysis',
            group: 'analysis_parameter',
            help: 'Setting this means Motif Discovery and Enrichment analysis will be performed.',
            example: true
        }
        insertsize: {
            description: 'Bowtie v1 maximum insert size (-X/--maxins <int>).',
            group: 'analysis_parameter',
            help: 'Specify maximum insert size for paired-end alignment (default: 600).',
            example: 600
        }
	    strandedness: {
            description: 'Bowtie v1 mate orientation (--fr/--rf/--ff).',
            group: 'analysis_parameter',
            help: 'The upstream/downstream mate orientation for paired-end alignment (default: --fr).',
            example: 'fr'
        }
    }

### ------------------------------------------------- ###
### ---------------- S E C T I O N 1 ---------------- ###
### ----------- Pre-process Analysis Files ---------- ###
### ------------------------------------------------- ###

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
            File R1end = select_first([fastqdump.R1end, string_sra[0]])
            File R2end = select_first([fastqdump.R2end, string_sra[0]])
        } # end scatter each sra

        Array[File] sample_R1_srafile = R1end
        Array[File] sample_R2_srafile = R2end
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
    #2. Make sure indexes are six else build indexes
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
    Array[File] actual_bowtie_index = select_first([bowtie_idx_2.bowtie_indexes, bowtie_idx.bowtie_indexes, bowtie_index])

    # FASTA faidx and chromsizes and effective genome size
    call samtools.faidx as samtools_faidx {
        # create FASTA index and chrom sizes files
        input :
            reference=reference
    }
    call util.effective_genome_size as egs {
        # effective genome size for FASTA
        input :
            reference=reference
    }

    # Process FASTQs
    if ( defined(sample_R1_fastq) ) {
        Array[String] string_fastq = [1] #buffer to allow for fastq optionality
        Array[File] s_R1_fastq = select_first([sample_R1_fastq, string_fastq])
        Array[File] sample_R1_fastqfile = s_R1_fastq
        Array[File] s_R2_fastq = select_first([sample_R2_fastq, string_fastq])
        Array[File] sample_R2_fastqfile = s_R2_fastq
        if (length(sample_R1_fastqfile) > 1) {
            call peaseq_util.sortfiles as R1_sorted { input: fastqfiles=sample_R1_fastqfile }
            call peaseq_util.sortfiles as R2_sorted { input: fastqfiles=sample_R2_fastqfile }
        }
        Array[File] sample_R1_fastqfiles = select_first([R1_sorted.allfiles, sample_R1_fastqfile])
        Array[File] sample_R2_fastqfiles = select_first([R2_sorted.allfiles, sample_R2_fastqfile])
    }

    # collate all fastqfiles
    Array[File] sample_R1 = flatten(select_all([sample_R1_srafile, sample_R1_fastqfiles]))
    Array[File] sample_R2 = flatten(select_all([sample_R2_srafile, sample_R2_fastqfiles]))
    Array[File] all_sample_fastqfiles = flatten(select_all([sample_R1_srafile, sample_R1_fastqfiles,sample_R2_srafile, sample_R2_fastqfiles]))

    # transpose to paired-end tuples
    Array[Pair[File, File]] sample_fastqfiles = zip(sample_R1, sample_R2)

    # if multiple fastqfiles are provided
    Boolean multi_fastqpair = if length(sample_fastqfiles) > 1 then true else false
    Boolean one_fastqpair = if length(sample_fastqfiles) == 1 then true else false

### ------------------------------------------------- ###
### ---------------- S E C T I O N 2 ---------------- ###
### ---- Single End (SE) Mode for all fastqfiles ---- ###
### ------------------------------------------------- ###

    scatter (eachfastq in all_sample_fastqfiles) {
        call fastqc.fastqc as indv_fastqc {
            input :
                inputfile=eachfastq,
                default_location=if multi_fastqpair then 'SAMPLE/' + sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/FastQC' else if defined(results_name) then results_name + '/single-end_mode/QC/FastQC' else sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/QC/FastQC'
        }

        call util.basicfastqstats as indv_bfs {
            input :
                fastqfile=eachfastq,
                default_location=if multi_fastqpair then 'SAMPLE/' + sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/SummaryStats' else if defined(results_name) then results_name + '/single-end_mode/QC/SummaryStats' else sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/QC/SummaryStats'
        }

        call mapping.mapping as indv_mapping {
            input :
                fastqfile=eachfastq,
                index_files=actual_bowtie_index,
                metricsfile=indv_bfs.metrics_out,
                blacklist=blacklist,
                default_location=if multi_fastqpair then 'SAMPLE/' + sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/BAM_files' else if defined(results_name) then results_name + '/single-end_mode/BAM_files' else sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/BAM_files'
        }

        call fastqc.fastqc as indv_bamfqc {
            input :
                inputfile=indv_mapping.sorted_bam,
                default_location=if multi_fastqpair then 'SAMPLE/' + sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/FastQC' else if defined(results_name) then results_name + '/single-end_mode/QC/FastQC' else sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/QC/FastQC'
        }

        call runspp.runspp as indv_runspp {
            input:
                bamfile=select_first([indv_mapping.bklist_bam, indv_mapping.sorted_bam])
        }

        call bedtools.bamtobed as indv_bamtobed {
            input:
                bamfile=select_first([indv_mapping.bklist_bam, indv_mapping.sorted_bam])
        }

        call util.evalstats as indv_summarystats {
            input:
                fastq_type="PEAseq Sample FASTQ",
                bambed=indv_bamtobed.bedfile,
                sppfile=indv_runspp.spp_out,
                fastqczip=indv_fastqc.zipfile,
                bamflag=indv_mapping.bam_stats,
                rmdupflag=indv_mapping.mkdup_stats,
                bkflag=indv_mapping.bklist_stats,
                fastqmetrics=indv_bfs.metrics_out,
                default_location=if multi_fastqpair then 'SAMPLE/' + sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/SummaryStats' else if defined(results_name) then results_name + '/single-end_mode/QC/SummaryStats' else sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/QC/SummaryStats'
        }
    } # end scatter (for eachfastq)

    # MERGE BAM files
    # Execute analysis on merge bam file
    # Analysis executed:
    #   Merge BAM (SE mode : for each fastq paired-end)

    call util.mergehtml {
        input:
            htmlfiles=indv_summarystats.xhtml,
            txtfiles=indv_summarystats.textfile,
            default_location='SAMPLE',
            outputfile = 'AllSamples-summary-stats.html'
    }

    call samtools.mergebam as SE_mergebam {
        input:
            bamfiles=indv_mapping.sorted_bam,
            metricsfiles=indv_bfs.metrics_out,
            default_location = if defined(results_name) then results_name + '/single-end_mode/BAM_files' else if multi_fastqpair then 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs' + '/single-end_mode/BAM_files' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/BAM_files',
            outputfile = 'AllMapped_' + length(all_sample_fastqfiles) + 'fastqs_SE.sorted.bam'
    }

    call fastqc.fastqc as SE_mergebamfqc {
        input:
            inputfile=SE_mergebam.mergebam,
            default_location=if defined(results_name) then results_name + '/single-end_mode/QC/FastQC' else if multi_fastqpair then 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs' + '/single-end_mode/QC/FastQC' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/QC/FastQC'
    }

    call samtools.indexstats as SE_mergeindexstats {
        input:
            bamfile=SE_mergebam.mergebam,
            default_location=if defined(results_name) then results_name + '/single-end_mode/BAM_files' else if multi_fastqpair then 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs' + '/single-end_mode/BAM_files' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/BAM_files'
    }

    if ( defined(blacklist) ) {
        # remove blacklist regions
        String string_blacklist = "" #buffer to allow for blacklist optionality
        File blacklist_file = select_first([blacklist, string_blacklist])
        call bedtools.intersect as SE_merge_rmblklist {
            input :
                fileA=SE_mergebam.mergebam,
                fileB=blacklist_file,
                default_location=if defined(results_name) then results_name + '/single-end_mode/BAM_files' else if multi_fastqpair then 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs' + '/single-end_mode/BAM_files' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/BAM_files',
                nooverlap=true
        }
        call samtools.indexstats as SE_merge_bklist {
            input :
                bamfile=SE_merge_rmblklist.intersect_out,
                default_location=if defined(results_name) then results_name + '/single-end_mode/BAM_files' else if multi_fastqpair then 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs' + '/single-end_mode/BAM_files' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/BAM_files'
        }
    } # end if blacklist provided

    File SE_mergebam_afterbklist = select_first([SE_merge_rmblklist.intersect_out, SE_mergebam.mergebam])

    call samtools.markdup as SE_merge_markdup {
        input :
            bamfile=SE_mergebam_afterbklist,
            default_location=if defined(results_name) then results_name + '/single-end_mode/BAM_files' else if multi_fastqpair then 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs' + '/single-end_mode/BAM_files' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/BAM_files'
    }

    call samtools.indexstats as SE_merge_mkdup {
        input :
            bamfile=SE_merge_markdup.mkdupbam,
            default_location=if defined(results_name) then results_name + '/single-end_mode/BAM_files' else if multi_fastqpair then 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs' + '/single-end_mode/BAM_files' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/BAM_files'
    }

### ------------------------------------------------- ###
### ---------------- S E C T I O N 3 ---------------- ###
### ---- Paired End (PE) Mode for all fastqfiles ---- ###
### ---- A: analysis if multiple FASTQs provided ---- ###
### ------------------------------------------------- ###

    if ( multi_fastqpair ) {
        scatter (fastqpair in sample_fastqfiles) {
            # Execute analysis on each fastq file provided
            # Analysis executed:
            #   Reference Alignment using Bowtie (-k2 -m2)
            #   Convert SAM to BAM
            #   FastQC on BAM files
            #   Remove Blacklists (if provided)
            #   Remove read duplicates
            #   Summary statistics on FASTQs
            #   Combine html files into one for easy viewing

            call util.basicfastqstats as indv_R1_bfs {
                input :
                    fastqfile=fastqpair.left,
                    default_location='SAMPLE/' + sub(basename(fastqpair.left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/SummaryStats'
            }
            call mapping.mapping as indv_PE_mapping {
                input :
                    fastqfile=fastqpair.left,
                    fastqfile_R2=fastqpair.right,
                    metricsfile=indv_R1_bfs.metrics_out,
                    insert_size=insertsize,
                    strandedness=strandedness,
                    index_files=actual_bowtie_index,
                    blacklist=blacklist,
                    paired_end=true,
                    default_location='SAMPLE/' + sub(basename(fastqpair.left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/BAM_files'
            }

            call fastqc.fastqc as indv_PE_bamfqc {
                input :
                    inputfile=indv_PE_mapping.sorted_bam,
                    default_location='SAMPLE/' + sub(basename(fastqpair.left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/FastQC'
            }

            call runspp.runspp as indv_PE_runspp {
                input:
                    bamfile=select_first([indv_PE_mapping.bklist_bam, indv_PE_mapping.sorted_bam])
            }

            call peaseq_util.pe_bamtobed as indv_PE_bamtobed {
                input:
                    bamfile=select_first([indv_PE_mapping.bklist_bam, indv_PE_mapping.sorted_bam])
            }

            call util.evalstats as indv_PE_summarystats {
                input:
                    fastq_type="PEAseq Sample FASTQ",
                    bambed=indv_PE_bamtobed.bedfile,
                    sppfile=indv_PE_runspp.spp_out,
                    fastqczip=indv_PE_bamfqc.zipfile,
                    bamflag=indv_PE_mapping.bam_stats,
                    rmdupflag=indv_PE_mapping.mkdup_stats,
                    bkflag=indv_PE_mapping.bklist_stats,
                    default_location='SAMPLE/' + sub(basename(fastqpair.left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/SummaryStats'
            }
        } # end scatter (for each sample fastq)

        # MERGE BAM FILES
        # Execute analysis on merge bam file
        # Analysis executed:
        #   Merge BAM (if more than 1 fastq is provided)
        #   FastQC on Merge BAM (AllMapped_<number>_mapped)

        # merge bam files and perform fasTQC if more than one is provided
        call util.mergehtml as PE_mergehtml {
            input:
                htmlfiles=indv_PE_summarystats.xhtml,
                txtfiles=indv_PE_summarystats.textfile,
                default_location='SAMPLE',
                outputfile = 'AllMapped_PEmode_peaseq-summary-stats.html'
        }

        call peaseq_util.pe_mergehtml as final_mergehtml {
            input:
                pe_htmlfiles=indv_PE_summarystats.xhtml,
                pe_txtfiles=indv_PE_summarystats.textfile,
                se_htmlfiles=indv_summarystats.xhtml,
                se_txtfiles=indv_summarystats.textfile,
                default_location='SAMPLE',
                outputfile = 'AllSamples-summary-stats.html'

        }

        call samtools.mergebam as PE_mergebam {
            input:
                bamfiles=indv_PE_mapping.as_sortedbam,
                metricsfiles=indv_bfs.metrics_out,
                paired_end=true,
                default_location = if defined(results_name) then results_name + '/BAM_files' else 'AllMapped_' + length(indv_PE_mapping.sorted_bam) + 'fastqpairs/BAM_files',
                outputfile = if defined(results_name) then results_name + '.sorted.bam' else 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs.sorted.bam',
                fixmatefile = if defined(results_name) then results_name + '.fixmate.bam' else 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs.fixmate.bam'
        }

        call fastqc.fastqc as PE_mergebamfqc {
            input:
	        inputfile=PE_mergebam.mergebam,
                default_location=sub(basename(PE_mergebam.mergebam),'.sorted.b.*$','') + '/QC/FastQC'
        }

        call samtools.indexstats as PE_mergeindexstats {
            input:
                bamfile=PE_mergebam.mergebam,
                default_location=sub(basename(PE_mergebam.mergebam),'.sorted.b.*$','') + '/BAM_files'
        }

        if ( defined(blacklist) ) {
            # remove blacklist regions
            String string_pe_blacklist = "" #buffer to allow for blacklist optionality
            File blacklist_pe = select_first([blacklist, string_pe_blacklist])
            String string_pe_fixmate = "" #buffer to allow for blacklist optionality
            File fixmate_pe = select_first([PE_mergebam.fixmatemergebam, string_pe_fixmate])
            call bedtools.pairtobed as PE_merge_rmblklist {
                input :
                    fileA=fixmate_pe,
                    fileB=blacklist_pe,
                    default_location=sub(basename(PE_mergebam.mergebam),'.sorted.b.*$','') + '/BAM_files'
            }
            call samtools.indexstats as PE_merge_bklist {
                input :
                    bamfile=PE_merge_rmblklist.pairtobed_out,
                    default_location=sub(basename(PE_mergebam.mergebam),'.sorted.b.*$','') + '/BAM_files'
            }
        } # end if blacklist provided

        File PE_mergebam_afterbklist = select_first([PE_merge_rmblklist.pairtobed_out, PE_mergebam.mergebam])

        call samtools.markdup as PE_merge_markdup {
            input :
                bamfile=PE_mergebam_afterbklist,
                default_location=sub(basename(PE_mergebam.mergebam),'.sorted.b.*$','') + '/BAM_files'
        }

        call samtools.indexstats as PE_merge_mkdup {
            input :
                bamfile=PE_merge_markdup.mkdupbam,
                default_location=sub(basename(PE_mergebam.mergebam),'.sorted.b.*$','') + '/BAM_files'
        }
    } # end if length(fastqfiles) > 1: multi_fastqpair

### ------------------------------------------------- ###
### ---------------- S E C T I O N 3 ---------------- ###
### ---- Paired End (PE) Mode for all fastqfiles ---- ###
### ------- B: analysis if one FASTQ provided ------- ###
### ------------------------------------------------- ###

    # if only one fastqfile is provided
    if ( one_fastqpair ) {
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

        call util.basicfastqstats as R1_bfs {
            input :
                fastqfile=sample_fastqfiles[0].left,
                default_location=if defined(results_name) then results_name + '/BAM_files' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/SummaryStats'
        }
        call mapping.mapping as uno_PE_mapping {
            input :
                fastqfile=sample_fastqfiles[0].left,
                fastqfile_R2=sample_fastqfiles[0].right,
                metricsfile=R1_bfs.metrics_out,
                insert_size=insertsize,
                strandedness=strandedness,
                index_files=actual_bowtie_index,
                blacklist=blacklist,
                paired_end=true,
                results_name=results_name,
                default_location=if defined(results_name) then results_name + '/BAM_files' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/BAM_files'
        }

        call fastqc.fastqc as uno_PE_bamfqc {
            input :
                inputfile=uno_PE_mapping.sorted_bam,
                default_location=sub(basename(uno_PE_mapping.sorted_bam),'.sorted.bam','') + '/QC/FastQC'
        }

        call runspp.runspp as uno_PE_runspp {
            input:
                bamfile=select_first([uno_PE_mapping.bklist_bam, uno_PE_mapping.sorted_bam])
        }

        call peaseq_util.pe_bamtobed as uno_PE_bamtobed {
            input:
                bamfile=select_first([uno_PE_mapping.bklist_bam, uno_PE_mapping.sorted_bam])
        }
    } # end if length(fastqfiles) == 1: one_fastqpair

### ------------------------------------------------- ###
### ---------------- S E C T I O N 4 ---------------- ###
### --------------- ChIP-seq Analysis --------------- ###
### ------------ A: analysis for SE mode  ----------- ###
### ------------------------------------------------- ###

    # ChIP-seq and downstream analysis
    # Execute analysis on merge bam file
    # Analysis executed:
    #   Peaks identification (SICER, MACS, ROSE)
    #   Motif analysis

    call macs.macs as SE_macs {
        input :
            bamfile=SE_mergebam_afterbklist,
            pvalue="1e-9",
            keep_dup="auto",
            egs=egs.genomesize,
            default_location = if defined(results_name) then results_name + '/single-end_mode/PEAKS/NARROW_peaks' + '/' + basename(SE_mergebam_afterbklist,'.bam') + '-p9_kd-auto' else if multi_fastqpair then 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs' + '/single-end_mode/PEAKS/NARROW_peaks' + '/' + basename(SE_mergebam_afterbklist,'.bam') + '-p9_kd-auto' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/PEAKS/NARROW_peaks' + '/' + basename(SE_mergebam_afterbklist,'.bam') + '-p9_kd-auto',
            coverage_location = if defined(results_name) then results_name + '/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + basename(SE_mergebam_afterbklist,'.bam') + '_p9_kd-auto' else if multi_fastqpair then 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs' + '/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + basename(SE_mergebam_afterbklist,'.bam') + '_p9_kd-auto' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + basename(SE_mergebam_afterbklist,'.bam') + '_p9_kd-auto'
    }

    call util.addreadme as SE_addreadme {
        input :
            default_location = if defined(results_name) then results_name + '/single-end_mode/PEAKS' else if multi_fastqpair then 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs' + '/single-end_mode/PEAKS' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/PEAKS'
    }

    call macs.macs as SE_all {
        input :
            bamfile=SE_mergebam_afterbklist,
            pvalue="1e-9",
            keep_dup="all",
            egs=egs.genomesize,
            default_location = if defined(results_name) then results_name + '/single-end_mode/PEAKS/NARROW_peaks' + '/' + basename(SE_mergebam_afterbklist,'.bam') + '-p9_kd-all' else if multi_fastqpair then 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs' + '/single-end_mode/PEAKS/NARROW_peaks' + '/' + basename(SE_mergebam_afterbklist,'.bam') + '-p9_kd-all' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/PEAKS/NARROW_peaks' + '/' + basename(SE_mergebam_afterbklist,'.bam') + '-p9_kd-all',
            coverage_location = if defined(results_name) then results_name + '/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + basename(SE_mergebam_afterbklist,'.bam') + '_p9_kd-all' else if multi_fastqpair then 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs' + '/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + basename(SE_mergebam_afterbklist,'.bam') + '_p9_kd-all' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + basename(SE_mergebam_afterbklist,'.bam') + '_p9_kd-all'
    }

    call macs.macs as SE_nomodel {
        input :
            bamfile=SE_mergebam_afterbklist,
            nomodel=true,
            egs=egs.genomesize,
            default_location = if defined(results_name) then results_name + '/single-end_mode/PEAKS/NARROW_peaks' + '/' + basename(SE_mergebam_afterbklist,'.bam') + '-nm' else if multi_fastqpair then 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs' + '/single-end_mode/PEAKS/NARROW_peaks' + '/' + basename(SE_mergebam_afterbklist,'.bam') + '-nm' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/PEAKS/NARROW_peaks' + '/' + basename(SE_mergebam_afterbklist,'.bam') + '-nm',
            coverage_location = if defined(results_name) then results_name + '/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + basename(SE_mergebam_afterbklist,'.bam') + '_nm' else if multi_fastqpair then 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs' + '/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + basename(SE_mergebam_afterbklist,'.bam') + '_nm' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + basename(SE_mergebam_afterbklist,'.bam') + '_nm'
    }

    call bamtogff.bamtogff as SE_bamtogff {
        input :
            gtffile=gtf,
            chromsizes=samtools_faidx.chromsizes,
            bamfile=SE_merge_markdup.mkdupbam,
            bamindex=SE_merge_mkdup.indexbam,
            default_location = if defined(results_name) then results_name + '/single-end_mode/BAM_Density' else if multi_fastqpair then 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs' + '/single-end_mode/BAM_Density' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/BAM_Density'
    }

    call bedtools.bamtobed as SE_forsicerbed {
        input :
            bamfile=SE_merge_markdup.mkdupbam
    }

    call sicer.sicer as SE_sicer {
        input :
            bedfile=SE_forsicerbed.bedfile,
            chromsizes=samtools_faidx.chromsizes,
            genome_fraction=egs.genomefraction,
            fragmentlength=SE_mergebam.avg_readlength,
            default_location = if defined(results_name) then results_name + '/single-end_mode/PEAKS/BROAD_peaks' else if multi_fastqpair then 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs' + '/single-end_mode/PEAKS/BROAD_peaks' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/PEAKS/BROAD_peaks',
            coverage_location = if defined(results_name) then results_name + '/single-end_mode/COVERAGE_files/BROAD_peaks' else if multi_fastqpair then 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs' + '/single-end_mode/COVERAGE_files/BROAD_peaks' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/COVERAGE_files/BROAD_peaks'
    }

    call rose.rose as SE_rose {
        input :
            gtffile=gtf,
            bamfile=SE_merge_markdup.mkdupbam,
            bamindex=SE_merge_mkdup.indexbam,
            bedfile_auto=SE_macs.peakbedfile,
            bedfile_all=SE_all.peakbedfile,
            default_location = if defined(results_name) then results_name + '/single-end_mode/PEAKS/STITCHED_peaks' else if multi_fastqpair then 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs' + '/single-end_mode/PEAKS/STITCHED_peaks' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/PEAKS/STITCHED_peaks'
    }

    call runspp.runspp as SE_runspp {
        input:
            bamfile=SE_mergebam_afterbklist
    }

    call util.peaksanno as SE_peaksanno {
        input :
            gtffile=gtf,
            bedfile=SE_macs.peakbedfile,
            chromsizes=samtools_faidx.chromsizes,
            summitfile=SE_macs.summitsfile,
            default_location = if defined(results_name) then results_name + '/single-end_mode/PEAKS_Annotation/NARROW_peaks' + '/' + sub(basename(SE_macs.peakbedfile),'_peaks.bed','') else if multi_fastqpair then 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs' + '/single-end_mode/PEAKS_Annotation/NARROW_peaks' + '/' + sub(basename(SE_macs.peakbedfile),'_peaks.bed','') else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/PEAKS_Annotation/NARROW_peaks' + '/' + sub(basename(SE_macs.peakbedfile),'_peaks.bed','')
    }

    call util.peaksanno as SE_all_peaksanno {
        input :
            gtffile=gtf,
            bedfile=SE_all.peakbedfile,
            chromsizes=samtools_faidx.chromsizes,
            summitfile=SE_all.summitsfile,
            default_location = if defined(results_name) then results_name + '/single-end_mode/PEAKS_Annotation/NARROW_peaks' + '/' + sub(basename(SE_all.peakbedfile),'_peaks.bed','') else if multi_fastqpair then 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs' + '/single-end_mode/PEAKS_Annotation/NARROW_peaks' + '/' + sub(basename(SE_all.peakbedfile),'_peaks.bed','') else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/PEAKS_Annotation/NARROW_peaks' + '/' + sub(basename(SE_all.peakbedfile),'_peaks.bed','')
    }

    call util.peaksanno as SE_nomodel_peaksanno {
        input :
            gtffile=gtf,
            bedfile=SE_nomodel.peakbedfile,
            chromsizes=samtools_faidx.chromsizes,
            summitfile=SE_nomodel.summitsfile,
            default_location = if defined(results_name) then results_name + '/single-end_mode/PEAKS_Annotation/NARROW_peaks' + '/' + sub(basename(SE_nomodel.peakbedfile),'_peaks.bed','') else if multi_fastqpair then 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs' + '/single-end_mode/PEAKS_Annotation/NARROW_peaks' + '/' + sub(basename(SE_nomodel.peakbedfile),'_peaks.bed','') else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/PEAKS_Annotation/NARROW_peaks' + '/' + sub(basename(SE_nomodel.peakbedfile),'_peaks.bed','')
    }

    call util.peaksanno as SE_sicer_peaksanno {
        input :
            gtffile=gtf,
            bedfile=SE_sicer.scoreisland,
            chromsizes=samtools_faidx.chromsizes,
            default_location = if defined(results_name) then results_name + '/single-end_mode/PEAKS_Annotation/BROAD_peaks' else if multi_fastqpair then 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs' + '/single-end_mode/PEAKS_Annotation/BROAD_peaks' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/PEAKS_Annotation/BROAD_peaks'
    }

    # Motif Analysis
    if (run_motifs) {
        call motifs.motifs as SE_motifs {
            input:
                reference=reference,
                reference_index=samtools_faidx.faidx_file,
                bedfile=SE_macs.peakbedfile,
                motif_databases=motif_databases,
                default_location = if defined(results_name) then results_name + '/single-end_mode/MOTIFS' else if multi_fastqpair then 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs' + '/single-end_mode/MOTIFS' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/MOTIFS'
        }

        call util.flankbed as SE_flankbed {
            input :
                bedfile=SE_macs.summitsfile,
                default_location = if defined(results_name) then results_name + '/single-end_mode/MOTIFS' else if multi_fastqpair then 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs' + '/single-end_mode/MOTIFS' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/MOTIFS'
        }

        call motifs.motifs as SE_flank {
            input:
                reference=reference,
                reference_index=samtools_faidx.faidx_file,
                bedfile=SE_flankbed.flankbedfile,
                motif_databases=motif_databases,
                default_location = if defined(results_name) then results_name + '/single-end_mode/MOTIFS' else if multi_fastqpair then 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs' + '/single-end_mode/MOTIFS' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/MOTIFS'
        }
    }

    call viz.visualization as SE_visualization {
        input:
            wigfile=SE_macs.wigfile,
            chromsizes=samtools_faidx.chromsizes,
            xlsfile=SE_macs.peakxlsfile,
            default_location = if defined(results_name) then results_name + '/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + sub(basename(SE_macs.peakbedfile),'_peaks.bed','') else if multi_fastqpair then 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs' + '/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + sub(basename(SE_macs.peakbedfile),'_peaks.bed','') else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + sub(basename(SE_macs.peakbedfile),'_peaks.bed','')
    }

    call viz.visualization as SE_vizall {
        input:
            wigfile=SE_all.wigfile,
            chromsizes=samtools_faidx.chromsizes,
            xlsfile=SE_all.peakxlsfile,
            default_location = if defined(results_name) then results_name + '/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + sub(basename(SE_all.peakbedfile),'_peaks.bed','') else if multi_fastqpair then 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs' + '/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + sub(basename(SE_all.peakbedfile),'_peaks.bed','') else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + sub(basename(SE_all.peakbedfile),'_peaks.bed','')
    }

    call viz.visualization as SE_viznomodel {
        input:
            wigfile=SE_nomodel.wigfile,
            chromsizes=samtools_faidx.chromsizes,
            xlsfile=SE_nomodel.peakxlsfile,
            default_location = if defined(results_name) then results_name + '/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + sub(basename(SE_nomodel.peakbedfile),'_peaks.bed','') else if multi_fastqpair then 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs' + '/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + sub(basename(SE_nomodel.peakbedfile),'_peaks.bed','') else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + sub(basename(SE_nomodel.peakbedfile),'_peaks.bed','')
    }

    call viz.visualization as SE_vizsicer {
        input:
            wigfile=SE_sicer.wigfile,
            chromsizes=samtools_faidx.chromsizes,
            default_location = if defined(results_name) then results_name + '/single-end_mode/COVERAGE_files/BROAD_peaks' else if multi_fastqpair then 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs' + '/single-end_mode/COVERAGE_files/BROAD_peaks' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/COVERAGE_files/BROAD_peaks'
    }

    call bedtools.bamtobed as SE_finalbed {
        input:
            bamfile=SE_mergebam_afterbklist
    }

    call sortbed.sortbed as SE_sortbed {
        input:
            bedfile=SE_finalbed.bedfile
    }

    call bedtools.intersect as SE_intersect {
        input:
            fileA=SE_macs.peakbedfile,
            fileB=SE_sortbed.sortbed_out,
            countoverlap=true,
            sorted=true
    }

### ------------------------------------------------- ###
### ---------------- S E C T I O N 4 ---------------- ###
### --------------- ChIP-seq Analysis --------------- ###
### ------------ B: analysis for PE mode  ----------- ###
### ------------------------------------------------- ###

    # ChIP-seq and downstream analysis
    # Execute analysis on merge bam file
    # Analysis executed:
    #   Peaks identification (SICER, MACS, ROSE)
    #   Motif analysis

    #collate correct files for downstream analysis
    File PE_sample_bam = select_first([PE_mergebam_afterbklist, uno_PE_mapping.bklist_bam, uno_PE_mapping.sorted_bam])

    call macs.macs as PE_macs {
        input :
            bamfile=PE_sample_bam,
            pvalue="1e-9",
            keep_dup="auto",
            egs=egs.genomesize,
            default_location=sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/PEAKS/NARROW_peaks' + '/' + basename(PE_sample_bam,'.bam') + '-p9_kd-auto',
            coverage_location=sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/COVERAGE_files/NARROW_peaks' + '/' + basename(PE_sample_bam,'.bam') + '_p9_kd-auto'
    }

    call util.addreadme as PE_addreadme {
        input :
            default_location=sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/PEAKS'
    }

    call macs.macs as PE_all {
        input :
            bamfile=PE_sample_bam,
            pvalue="1e-9",
            keep_dup="all",
            egs=egs.genomesize,
            default_location=sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/PEAKS/NARROW_peaks' + '/' + basename(PE_sample_bam,'.bam') + '-p9_kd-all',
            coverage_location=sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/COVERAGE_files/NARROW_peaks' + '/' + basename(PE_sample_bam,'.bam') + '_p9_kd-all'
    }

    call macs.macs as PE_nomodel {
        input :
            bamfile=PE_sample_bam,
            nomodel=true,
            egs=egs.genomesize,
            default_location=sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/PEAKS/NARROW_peaks' + '/' + basename(PE_sample_bam,'.bam') + '-nm',
            coverage_location=sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/COVERAGE_files/NARROW_peaks' + '/' + basename(PE_sample_bam,'.bam') + '_nm'
    }

    call peaseq_util.fraggraph {
        input :
            bamfile=select_first([PE_merge_markdup.mkdupbam, uno_PE_mapping.mkdup_bam]),
            chromsizes=samtools_faidx.chromsizes,
            default_location=sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/COVERAGE_files/BAM_Fragments',
            bam_location=sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/BAM_files',
            annotation_location=sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/QC/Fragments'
    }

    call samtools.indexstats as frag_index {
        input :
            bamfile=fraggraph.fragbamfile,
            default_location=sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/BAM_files'
    }

    call bamtogff.bamtogff as PE_bamtogff {
        input :
            gtffile=gtf,
            chromsizes=samtools_faidx.chromsizes,
            bamfile=fraggraph.fragbamfile,
            bamindex=frag_index.indexbam,
            default_location=sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/BAM_Density'
    }

    call sicer.sicer as PE_sicer {
        input :
            bedfile=fraggraph.bedpefile,
            paired_end=true,
            gap_size=600,
            chromsizes=samtools_faidx.chromsizes,
            genome_fraction=egs.genomefraction,
            default_location=sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/PEAKS/BROAD_peaks',
            coverage_location=sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/COVERAGE_files/BROAD_peaks'
    }

    call rose.rose as PE_rose {
        input :
            gtffile=gtf,
            bamfile=PE_sample_bam,
            bamindex=select_first([PE_merge_bklist.indexbam, PE_mergeindexstats.indexbam, uno_PE_mapping.bklist_index, uno_PE_mapping.bam_index]),
            bedfile_auto=PE_macs.peakbedfile,
            bedfile_all=PE_all.peakbedfile,
            default_location=sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/PEAKS/STITCHED_peaks'
    }

    call runspp.runspp as PE_runspp {
        input:
            bamfile=PE_sample_bam
    }

    call util.peaksanno as PE_peaksanno {
        input :
            gtffile=gtf,
            bedfile=PE_macs.peakbedfile,
            chromsizes=samtools_faidx.chromsizes,
            summitfile=PE_macs.summitsfile,
            default_location=sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/PEAKS_Annotation/NARROW_peaks' + '/' + sub(basename(PE_macs.peakbedfile),'_peaks.bed','')
    }

    call util.peaksanno as PE_all_peaksanno {
        input :
            gtffile=gtf,
            bedfile=PE_all.peakbedfile,
            chromsizes=samtools_faidx.chromsizes,
            summitfile=PE_all.summitsfile,
            default_location=sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/PEAKS_Annotation/NARROW_peaks' + '/' + sub(basename(PE_all.peakbedfile),'_peaks.bed','')
    }

    call util.peaksanno as PE_nomodel_peaksanno {
        input :
            gtffile=gtf,
            bedfile=PE_nomodel.peakbedfile,
            chromsizes=samtools_faidx.chromsizes,
            summitfile=PE_nomodel.summitsfile,
            default_location=sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/PEAKS_Annotation/NARROW_peaks' + '/' + sub(basename(PE_nomodel.peakbedfile),'_peaks.bed','')
    }

    call util.peaksanno as PE_sicer_peaksanno {
        input :
            gtffile=gtf,
            bedfile=PE_sicer.scoreisland,
            chromsizes=samtools_faidx.chromsizes,
            default_location=sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/PEAKS_Annotation/BROAD_peaks'
    }

    # Motif Analysis
    if (run_motifs) {
        call motifs.motifs as PE_motifs {
            input:
                reference=reference,
                reference_index=samtools_faidx.faidx_file,
                bedfile=PE_macs.peakbedfile,
                motif_databases=motif_databases,
                default_location=sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/MOTIFS'
        }

        call util.flankbed as PE_flankbed {
            input :
                bedfile=PE_macs.summitsfile,
                default_location=sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/MOTIFS'
        }

        call motifs.motifs as PE_flank {
            input:
                reference=reference,
                reference_index=samtools_faidx.faidx_file,
                bedfile=PE_flankbed.flankbedfile,
                motif_databases=motif_databases,
                default_location=sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/MOTIFS'
        }
    }

    call viz.visualization as PE_visualization {
        input:
            wigfile=PE_macs.wigfile,
            chromsizes=samtools_faidx.chromsizes,
            xlsfile=PE_macs.peakxlsfile,
            default_location=sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/COVERAGE_files/NARROW_peaks' + '/' + sub(basename(PE_macs.peakbedfile),'_peaks.bed','')
    }

    call viz.visualization as PE_vizall {
        input:
            wigfile=PE_all.wigfile,
            chromsizes=samtools_faidx.chromsizes,
            xlsfile=PE_all.peakxlsfile,
            default_location=sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/COVERAGE_files/NARROW_peaks' + '/' + sub(basename(PE_all.peakbedfile),'_peaks.bed','')
    }

    call viz.visualization as PE_viznomodel {
        input:
            wigfile=PE_nomodel.wigfile,
            chromsizes=samtools_faidx.chromsizes,
            xlsfile=PE_nomodel.peakxlsfile,
            default_location=sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/COVERAGE_files/NARROW_peaks' + '/' + sub(basename(PE_nomodel.peakbedfile),'_peaks.bed','')
    }

    call viz.visualization as PE_vizsicer {
        input:
            wigfile=PE_sicer.wigfile,
            chromsizes=samtools_faidx.chromsizes,
            default_location=sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/COVERAGE_files/BROAD_peaks'
    }

    call bedtools.bamtobed as PE_finalbed {
        input:
            bamfile=PE_sample_bam
    }

    call sortbed.sortbed as PE_sortbed {
        input:
            bedfile=PE_finalbed.bedfile
    }

    call bedtools.intersect as PE_intersect {
        input:
            fileA=PE_macs.peakbedfile,
            fileB=PE_sortbed.sortbed_out,
            countoverlap=true,
            sorted=true
    }

    # Final MERGE output file
    File mergehtmlfile =  select_first([final_mergehtml.mergefile, mergehtml.mergefile])
    File mergetxtfile =  select_first([final_mergehtml.mergetxt, mergehtml.mergetxt])

### ------------------------------------------------- ###
### ----------------- S E C T I O N 5 --------------- ###
### --------------- Summary Statistics -------------- ###
### ------------ A: analysis for SE mode  ----------- ###
### ------------------------------------------------- ###
    call util.evalstats as SE_summarystats {
        input:
            fastq_type="PEAseq SEmode Comprehensive",
            bambed=SE_finalbed.bedfile,
            sppfile=SE_runspp.spp_out,
            fastqczip=SE_mergebamfqc.zipfile,
            bamflag=SE_mergeindexstats.flagstats,
            rmdupflag=SE_merge_mkdup.flagstats,
            bkflag=SE_merge_bklist.flagstats,
            countsfile=SE_intersect.intersect_out,
            peaksxls=SE_macs.peakxlsfile,
            enhancers=SE_rose.enhancers,
            superenhancers=SE_rose.super_enhancers,
            default_location = if defined(results_name) then results_name + '/single-end_mode/QC/SummaryStats' else if multi_fastqpair then 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs' + '/single-end_mode/QC/SummaryStats' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/QC/SummaryStats'
    }

    call util.summaryreport as merge_overallsummary {
        input:
            sampleqc_html=mergehtml.xhtml,
            overallqc_html=SE_summarystats.xhtml,
            sampleqc_txt=mergehtml.mergetxt,
            overallqc_txt=SE_summarystats.textfile,
            default_location=if defined(results_name) then results_name + '/single-end_mode' else if multi_fastqpair then 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs' + '/single-end_mode' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode'
    }
### ------------------------------------------------- ###
### ----------------- S E C T I O N 5 --------------- ###
### --------------- Summary Statistics -------------- ###
### ------------ B: analysis for PE mode  ----------- ###
### ------------------------------------------------- ###

    String string_qual = "" #buffer to allow for optionality in if statement

    if ( one_fastqpair ) {
        call util.evalstats as PE_uno_summarystats {
            # SUMMARY STATISTICS of sample file (only 1 sample file provided)
            input:
                fastq_type="PEAseq PEmode Comprehensive",
                bambed=PE_finalbed.bedfile,
                sppfile=PE_runspp.spp_out,
                fastqczip=select_first([uno_PE_bamfqc.zipfile, string_qual]),
                bamflag=uno_PE_mapping.bam_stats,
                rmdupflag=uno_PE_mapping.mkdup_stats,
                bkflag=uno_PE_mapping.bklist_stats,
                fastqmetrics=indv_bfs.metrics_out[0],
                countsfile=PE_intersect.intersect_out,
                peaksxls=PE_macs.peakxlsfile,
                enhancers=PE_rose.enhancers,
                superenhancers=PE_rose.super_enhancers,
                default_location=sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/QC/SummaryStats'
        }

        call peaseq_util.pairedend_summaryreport as PE_uno_overallsummary {
            # Presenting all quality stats for the analysis
            input:
                sampleqc_se_html=mergehtml.xhtml,
                overallqc_se_html=SE_summarystats.xhtml,
                sampleqc_se_txt=mergehtml.mergetxt,
                overallqc_se_txt=SE_summarystats.textfile,
                overallqc_pe_html=PE_uno_summarystats.xhtml,
                overallqc_pe_txt=PE_uno_summarystats.textfile,
                outputfile = if defined(results_name) then results_name + '.peaseq_report.html' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '.peaseq_report.html'
        }
    } # end if one_fastqpair

    if ( multi_fastqpair ) {
        call util.evalstats as PE_merge_summarystats {
            # SUMMARY STATISTICS of all samples files (more than 1 sample file provided)
            input:
                fastq_type="PEAseq PEmode Comprehensive",
                bambed=PE_finalbed.bedfile,
                sppfile=PE_runspp.spp_out,
                fastqczip=select_first([PE_mergebamfqc.zipfile, string_qual]),
                bamflag=PE_mergeindexstats.flagstats,
                rmdupflag=PE_merge_mkdup.flagstats,
                bkflag=PE_merge_bklist.flagstats,
                countsfile=PE_intersect.intersect_out,
                peaksxls=PE_macs.peakxlsfile,
                enhancers=PE_rose.enhancers,
                superenhancers=PE_rose.super_enhancers,
                default_location=sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/QC/SummaryStats'
        }

        call peaseq_util.pairedend_summaryreport as PE_merge_overallsummary {
            # Presenting all quality stats for the analysis
            input:
                sampleqc_se_html=mergehtml.xhtml,
                overallqc_se_html=SE_summarystats.xhtml,
                sampleqc_se_txt=mergehtml.mergetxt,
                overallqc_se_txt=SE_summarystats.textfile,
                sampleqc_pe_html=PE_mergehtml.xhtml,
                overallqc_pe_html=PE_merge_summarystats.xhtml,
                sampleqc_pe_txt=PE_mergehtml.mergetxt,
                overallqc_pe_txt=PE_merge_summarystats.textfile,
                outputfile = if defined(results_name) then results_name + '.peaseq_report.html' else 'AllMapped_' + length(sample_fastqfiles) + 'fastqpairs.peaseq_report.html'
        }
    } # end if multi_fastqpair


### ------------------------------------------------- ###
### ---------------- S E C T I O N 6 ---------------- ###
### ------------------ OUTPUT FILES ----------------- ###
### ------------------------------------------------- ###

    output {

        #FASTQC
        Array[File?]? indv_s_htmlfile = indv_fastqc.htmlfile
        Array[File?]? indv_s_zipfile = indv_fastqc.zipfile
        Array[File?]? indv_s_bam_htmlfile = indv_bamfqc.htmlfile
        Array[File?]? indv_s_bam_zipfile = indv_bamfqc.zipfile
        File? s_mergebam_htmlfile = SE_mergebamfqc.htmlfile
        File? s_mergebam_zipfile = SE_mergebamfqc.zipfile

        Array[File?]? indv_sp_bam_htmlfile = indv_PE_bamfqc.htmlfile
        Array[File?]? indv_sp_bam_zipfile = indv_PE_bamfqc.zipfile
        File? sp_mergebam_htmlfile = PE_mergebamfqc.htmlfile
        File? sp_mergebam_zipfile = PE_mergebamfqc.zipfile

        File? uno_s_bam_htmlfile = uno_PE_bamfqc.htmlfile
        File? uno_s_bam_zipfile = uno_PE_bamfqc.zipfile

        #BASICMETRICS
        Array[File?]? s_metrics_out = indv_bfs.metrics_out

        #BAMFILES
        Array[File?]? indv_s_sortedbam = indv_mapping.sorted_bam
        Array[File?]? indv_s_indexbam = indv_mapping.bam_index
        Array[File?]? indv_s_bkbam = indv_mapping.bklist_bam
        Array[File?]? indv_s_bkindexbam = indv_mapping.bklist_index
        Array[File?]? indv_s_rmbam = indv_mapping.mkdup_bam
        Array[File?]? indv_s_rmindexbam = indv_mapping.mkdup_index
        Array[File?]? indv_sp_sortedbam = indv_PE_mapping.sorted_bam
        Array[File?]? indv_sp_indexbam = indv_PE_mapping.bam_index
        Array[File?]? indv_sp_bkbam = indv_PE_mapping.bklist_bam
        Array[File?]? indv_sp_bkindexbam = indv_PE_mapping.bklist_index
        Array[File?]? indv_sp_rmbam = indv_PE_mapping.mkdup_bam
        Array[File?]? indv_sp_rmindexbam = indv_PE_mapping.mkdup_index

        File? uno_s_sortedbam = uno_PE_mapping.sorted_bam
        File? uno_s_indexstatsbam = uno_PE_mapping.bam_index
        File? uno_s_bkbam = uno_PE_mapping.bklist_bam
        File? uno_s_bkindexbam = uno_PE_mapping.bklist_index
        File? uno_s_rmbam = uno_PE_mapping.mkdup_bam
        File? uno_s_rmindexbam = uno_PE_mapping.mkdup_index
        File? s_mergebamfile = SE_mergebam.mergebam
        File? s_mergebamindex = SE_mergeindexstats.indexbam
        File? s_bkbam = SE_merge_rmblklist.intersect_out
        File? s_bkindexbam = SE_merge_bklist.indexbam
        File? s_rmbam = SE_merge_markdup.mkdupbam
        File? s_rmindexbam = SE_merge_mkdup.indexbam

        File? sp_mergebamfile = PE_mergebam.mergebam
        File? sp_mergebamindex = PE_mergeindexstats.indexbam
        File? sp_bkbam = PE_merge_rmblklist.pairtobed_out
        File? sp_bkindexbam = PE_merge_bklist.indexbam
        File? sp_rmbam = PE_merge_markdup.mkdupbam
        File? sp_rmindexbam = PE_merge_mkdup.indexbam

        File? s_fragments_bam = fraggraph.fragbamfile
        File? s_fragments_indexbam = frag_index.indexbam

        #MACS
        File? peakbedfile = SE_macs.peakbedfile
        File? peakxlsfile = SE_macs.peakxlsfile
        File? summitsfile = SE_macs.summitsfile
        File? negativexlsfile = SE_macs.negativepeaks
        File? wigfile = SE_macs.wigfile
        File? all_peakbedfile = SE_all.peakbedfile
        File? all_peakxlsfile = SE_all.peakxlsfile
        File? all_summitsfile = SE_all.summitsfile
        File? all_negativexlsfile = SE_all.negativepeaks
        File? all_wigfile = SE_all.wigfile
        File? nm_peakbedfile = SE_nomodel.peakbedfile
        File? nm_peakxlsfile = SE_nomodel.peakxlsfile
        File? nm_summitsfile = SE_nomodel.summitsfile
        File? nm_negativexlsfile = SE_nomodel.negativepeaks
        File? nm_wigfile = SE_nomodel.wigfile
        File? readme_peaks = SE_addreadme.readme_peaks
        File? sp_peakbedfile = PE_macs.peakbedfile
        File? sp_peakxlsfile = PE_macs.peakxlsfile
        File? sp_summitsfile = PE_macs.summitsfile
        File? sp_negativexlsfile = PE_macs.negativepeaks
        File? sp_wigfile = PE_macs.wigfile
        File? sp_all_peakbedfile = PE_all.peakbedfile
        File? sp_all_peakxlsfile = PE_all.peakxlsfile
        File? sp_all_summitsfile = PE_all.summitsfile
        File? sp_all_negativexlsfile = PE_all.negativepeaks
        File? sp_all_wigfile = PE_all.wigfile
        File? sp_nm_peakbedfile = PE_nomodel.peakbedfile
        File? sp_nm_peakxlsfile = PE_nomodel.peakxlsfile
        File? sp_nm_summitsfile = PE_nomodel.summitsfile
        File? sp_nm_negativexlsfile = PE_nomodel.negativepeaks
        File? sp_nm_wigfile = PE_nomodel.wigfile
        File? sp_readme_peaks = PE_addreadme.readme_peaks

        #SICER
        File? scoreisland = SE_sicer.scoreisland
        File? sicer_wigfile = SE_sicer.wigfile
        File? sp_scoreisland = PE_sicer.scoreisland
        File? sp_sicer_wigfile = PE_sicer.wigfile

        #ROSE
        File? pngfile = SE_rose.pngfile
        File? mapped_union = SE_rose.mapped_union
        File? mapped_stitch = SE_rose.mapped_stitch
        File? enhancers = SE_rose.enhancers
        File? super_enhancers = SE_rose.super_enhancers
        File? gff_file = SE_rose.gff_file
        File? gff_union = SE_rose.gff_union
        File? union_enhancers = SE_rose.union_enhancers
        File? stitch_enhancers = SE_rose.stitch_enhancers
        File? e_to_g_enhancers = SE_rose.e_to_g_enhancers
        File? g_to_e_enhancers = SE_rose.g_to_e_enhancers
        File? e_to_g_super_enhancers = SE_rose.e_to_g_super_enhancers
        File? g_to_e_super_enhancers = SE_rose.g_to_e_super_enhancers
        File? sp_pngfile = PE_rose.pngfile
        File? sp_mapped_union = PE_rose.mapped_union
        File? sp_mapped_stitch = PE_rose.mapped_stitch
        File? sp_enhancers = PE_rose.enhancers
        File? sp_super_enhancers = PE_rose.super_enhancers
        File? sp_gff_file = PE_rose.gff_file
        File? sp_gff_union = PE_rose.gff_union
        File? sp_union_enhancers = PE_rose.union_enhancers
        File? sp_stitch_enhancers = PE_rose.stitch_enhancers
        File? sp_e_to_g_enhancers = PE_rose.e_to_g_enhancers
        File? sp_g_to_e_enhancers = PE_rose.g_to_e_enhancers
        File? sp_e_to_g_super_enhancers = PE_rose.e_to_g_super_enhancers
        File? sp_g_to_e_super_enhancers = PE_rose.g_to_e_super_enhancers

        #MOTIFS
        File? flankbedfile = SE_flankbed.flankbedfile
        File? ame_tsv = SE_motifs.ame_tsv
        File? ame_html = SE_motifs.ame_html
        File? ame_seq = SE_motifs.ame_seq
        File? meme = SE_motifs.meme_out
        File? meme_summary = SE_motifs.meme_summary
        File? summit_ame_tsv = SE_flank.ame_tsv
        File? summit_ame_html = SE_flank.ame_html
        File? summit_ame_seq = SE_flank.ame_seq
        File? summit_meme = SE_flank.meme_out
        File? summit_meme_summary = SE_flank.meme_summary
        File? sp_flankbedfile = PE_flankbed.flankbedfile
        File? sp_ame_tsv = PE_motifs.ame_tsv
        File? sp_ame_html = PE_motifs.ame_html
        File? sp_ame_seq = PE_motifs.ame_seq
        File? sp_meme = PE_motifs.meme_out
        File? sp_meme_summary = PE_motifs.meme_summary
        File? sp_summit_ame_tsv = PE_flank.ame_tsv
        File? sp_summit_ame_html = PE_flank.ame_html
        File? sp_summit_ame_seq = PE_flank.ame_seq
        File? sp_summit_meme = PE_flank.meme_out
        File? sp_summit_meme_summary = PE_flank.meme_summary

        #BAM2GFF
        File? s_matrices = SE_bamtogff.s_matrices
        File? densityplot = SE_bamtogff.densityplot
        File? pdf_gene = SE_bamtogff.pdf_gene
        File? pdf_h_gene = SE_bamtogff.pdf_h_gene
        File? png_h_gene = SE_bamtogff.png_h_gene
        File? jpg_h_gene = SE_bamtogff.jpg_h_gene
        File? pdf_promoters = SE_bamtogff.pdf_promoters
        File? pdf_h_promoters = SE_bamtogff.pdf_h_promoters
        File? png_h_promoters = SE_bamtogff.png_h_promoters
        File? jpg_h_promoters = SE_bamtogff.jpg_h_promoters

        File? sp_s_matrices = PE_bamtogff.s_matrices
        File? sp_densityplot = PE_bamtogff.densityplot
        File? sp_pdf_gene = PE_bamtogff.pdf_gene
        File? sp_pdf_h_gene = PE_bamtogff.pdf_h_gene
        File? sp_png_h_gene = PE_bamtogff.png_h_gene
        File? sp_jpg_h_gene = PE_bamtogff.jpg_h_gene
        File? sp_pdf_promoters = PE_bamtogff.pdf_promoters
        File? sp_pdf_h_promoters = PE_bamtogff.pdf_h_promoters
        File? sp_png_h_promoters = PE_bamtogff.png_h_promoters
        File? sp_jpg_h_promoters = PE_bamtogff.jpg_h_promoters

        #PEAKS-ANNOTATION
        File? peak_promoters = SE_peaksanno.peak_promoters
        File? peak_genebody = SE_peaksanno.peak_genebody
        File? peak_window = SE_peaksanno.peak_window
        File? peak_closest = SE_peaksanno.peak_closest
        File? peak_comparison = SE_peaksanno.peak_comparison
        File? gene_comparison = SE_peaksanno.gene_comparison
        File? pdf_comparison = SE_peaksanno.pdf_comparison
        File? all_peak_promoters = SE_all_peaksanno.peak_promoters
        File? all_peak_genebody = SE_all_peaksanno.peak_genebody
        File? all_peak_window = SE_all_peaksanno.peak_window
        File? all_peak_closest = SE_all_peaksanno.peak_closest
        File? all_peak_comparison = SE_all_peaksanno.peak_comparison
        File? all_gene_comparison = SE_all_peaksanno.gene_comparison
        File? all_pdf_comparison = SE_all_peaksanno.pdf_comparison
        File? nomodel_peak_promoters = SE_nomodel_peaksanno.peak_promoters
        File? nomodel_peak_genebody = SE_nomodel_peaksanno.peak_genebody
        File? nomodel_peak_window = SE_nomodel_peaksanno.peak_window
        File? nomodel_peak_closest = SE_nomodel_peaksanno.peak_closest
        File? nomodel_peak_comparison = SE_nomodel_peaksanno.peak_comparison
        File? nomodel_gene_comparison = SE_nomodel_peaksanno.gene_comparison
        File? nomodel_pdf_comparison = SE_nomodel_peaksanno.pdf_comparison
        File? sicer_peak_promoters = SE_sicer_peaksanno.peak_promoters
        File? sicer_peak_genebody = SE_sicer_peaksanno.peak_genebody
        File? sicer_peak_window = SE_sicer_peaksanno.peak_window
        File? sicer_peak_closest = SE_sicer_peaksanno.peak_closest
        File? sicer_peak_comparison = SE_sicer_peaksanno.peak_comparison
        File? sicer_gene_comparison = SE_sicer_peaksanno.gene_comparison
        File? sicer_pdf_comparison = SE_sicer_peaksanno.pdf_comparison

        File? sp_peak_promoters = PE_peaksanno.peak_promoters
        File? sp_peak_genebody = PE_peaksanno.peak_genebody
        File? sp_peak_window = PE_peaksanno.peak_window
        File? sp_peak_closest = PE_peaksanno.peak_closest
        File? sp_peak_comparison = PE_peaksanno.peak_comparison
        File? sp_gene_comparison = PE_peaksanno.gene_comparison
        File? sp_pdf_comparison = PE_peaksanno.pdf_comparison
        File? sp_all_peak_promoters = PE_all_peaksanno.peak_promoters
        File? sp_all_peak_genebody = PE_all_peaksanno.peak_genebody
        File? sp_all_peak_window = PE_all_peaksanno.peak_window
        File? sp_all_peak_closest = PE_all_peaksanno.peak_closest
        File? sp_all_peak_comparison = PE_all_peaksanno.peak_comparison
        File? sp_all_gene_comparison = PE_all_peaksanno.gene_comparison
        File? sp_all_pdf_comparison = PE_all_peaksanno.pdf_comparison
        File? sp_nomodel_peak_promoters = PE_nomodel_peaksanno.peak_promoters
        File? sp_nomodel_peak_genebody = PE_nomodel_peaksanno.peak_genebody
        File? sp_nomodel_peak_window = PE_nomodel_peaksanno.peak_window
        File? sp_nomodel_peak_closest = PE_nomodel_peaksanno.peak_closest
        File? sp_nomodel_peak_comparison = PE_nomodel_peaksanno.peak_comparison
        File? sp_nomodel_gene_comparison = PE_nomodel_peaksanno.gene_comparison
        File? sp_nomodel_pdf_comparison = PE_nomodel_peaksanno.pdf_comparison
        File? sp_sicer_peak_promoters = PE_sicer_peaksanno.peak_promoters
        File? sp_sicer_peak_genebody = PE_sicer_peaksanno.peak_genebody
        File? sp_sicer_peak_window = PE_sicer_peaksanno.peak_window
        File? sp_sicer_peak_closest = PE_sicer_peaksanno.peak_closest
        File? sp_sicer_peak_comparison = PE_sicer_peaksanno.peak_comparison
        File? sp_sicer_gene_comparison = PE_sicer_peaksanno.gene_comparison
        File? sp_sicer_pdf_comparison = PE_sicer_peaksanno.pdf_comparison

        #VISUALIZATION
        File? bigwig = SE_visualization.bigwig
        File? norm_wig = SE_visualization.norm_wig
        File? tdffile = SE_visualization.tdffile
        File? n_bigwig = SE_viznomodel.bigwig
        File? n_norm_wig = SE_viznomodel.norm_wig
        File? n_tdffile = SE_viznomodel.tdffile
        File? a_bigwig = SE_vizall.bigwig
        File? a_norm_wig = SE_vizall.norm_wig
        File? a_tdffile = SE_vizall.tdffile
        File? s_bigwig = SE_vizsicer.bigwig
        File? s_norm_wig = SE_vizsicer.norm_wig
        File? s_tdffile = SE_vizsicer.tdffile

        File? sp_bigwig = PE_visualization.bigwig
        File? sp_norm_wig = PE_visualization.norm_wig
        File? sp_tdffile = PE_visualization.tdffile
        File? sp_n_bigwig = PE_viznomodel.bigwig
        File? sp_n_norm_wig = PE_viznomodel.norm_wig
        File? sp_n_tdffile = PE_viznomodel.tdffile
        File? sp_a_bigwig = PE_vizall.bigwig
        File? sp_a_norm_wig = PE_vizall.norm_wig
        File? sp_a_tdffile = PE_vizall.tdffile
        File? sp_s_bigwig = PE_vizsicer.bigwig
        File? sp_s_norm_wig = PE_vizsicer.norm_wig
        File? sp_s_tdffile = PE_vizsicer.tdffile

        File? sf_bigwig = fraggraph.bigwigfile
        File? sf_tdffile = fraggraph.tdffile
        File? sf_wigfile = fraggraph.wigfile

        #QC-STATS
        Array[File?]? s_qc_statsfile = indv_summarystats.statsfile
        Array[File?]? s_qc_htmlfile = indv_summarystats.htmlfile
        Array[File?]? s_qc_textfile = indv_summarystats.textfile
        File? statsfile = SE_summarystats.statsfile
        File? htmlfile = SE_summarystats.htmlfile
        File? textfile = SE_summarystats.textfile

        File? s_summaryhtml = merge_overallsummary.summaryhtml
        File? s_summarytxt = merge_overallsummary.summarytxt

        File? s_qc_mergehtml = final_mergehtml.mergefile

        Array[File?]? sp_qc_statsfile = indv_PE_summarystats.statsfile
        Array[File?]? sp_qc_htmlfile = indv_PE_summarystats.htmlfile
        Array[File?]? sp_qc_textfile = indv_PE_summarystats.textfile
        File? s_uno_statsfile = PE_uno_summarystats.statsfile
        File? s_uno_htmlfile = PE_uno_summarystats.htmlfile
        File? s_uno_textfile = PE_uno_summarystats.textfile
        File? sp_statsfile = PE_merge_summarystats.statsfile
        File? sp_htmlfile = PE_merge_summarystats.htmlfile
        File? sp_textfile = PE_merge_summarystats.textfile
        File? summaryhtml = PE_uno_overallsummary.summaryhtml
        File? summarytxt = PE_uno_overallsummary.summarytxt
        File? sp_summaryhtml = PE_merge_overallsummary.summaryhtml
        File? sp_summarytxt = PE_merge_overallsummary.summarytxt
        File? s_fragsize = fraggraph.fragsizepng
    }
}
