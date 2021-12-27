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

workflow peaseq {
    String pipeline_ver = 'v2.0.0'

    meta {
        title: 'PEAseq Analysis'
        summary: 'Paired-End Antibody Sequencing (PEAseq) Pipeline'
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
                    changes: "Initial release"
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
        Int? bowtie_insert_size = 600
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
        bowtie_insert_size: {
            description: 'Bowtie v1 maximum insert size (-X/--maxins <int>).',
            group: 'analysis_parameter',
            help: 'Specify maximum insert size for paired-end alignment (default: 600).',
            example: 600
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

        Array[File] sample_R1_srafile_ = R1end 
        Array[File] sample_R2_srafile_ = R2end 
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
    Array[File] bowtie_index_ = select_first([bowtie_idx_2.bowtie_indexes, bowtie_idx.bowtie_indexes, bowtie_index])

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
        Array[File] sample_R1_fastqfile_ = s_R1_fastq
        Array[File] s_R2_fastq = select_first([sample_R2_fastq, string_fastq])
        Array[File] sample_R2_fastqfile_ = s_R2_fastq
    } 

    # collate all fastqfiles
    Array[File] sample_R1 = flatten(select_all([sample_R1_srafile_, sample_R1_fastqfile_]))
    Array[File] sample_R2 = flatten(select_all([sample_R2_srafile_, sample_R2_fastqfile_]))
    Array[File] all_sample_fastqfiles = flatten(select_all([sample_R1_srafile_, sample_R1_fastqfile_,sample_R2_srafile_, sample_R2_fastqfile_]))

    # transpose to paired-end tuples
    Array[Pair[File, File]] sample_fastqfiles = zip(sample_R1, sample_R2)

### ------------------------------------------------- ###
### ---------------- S E C T I O N 2 ---------------- ###
### ---- Single End (SE) Mode for all fastqfiles ---- ###
### ------------------------------------------------- ###

    scatter (eachfastq in all_sample_fastqfiles) {
        call fastqc.fastqc as indv_fastqc {
            input :
                inputfile=eachfastq,
                default_location='SAMPLE/' + sub(basename(eachfastq),'_R?[12].*\.f.*q\.gz','') + '/QC/FastQC'
        }

        call util.basicfastqstats as indv_bfs {
            input :
                fastqfile=eachfastq,
                default_location='SAMPLE/' + sub(basename(eachfastq),'_R?[12].*\.f.*q\.gz','') + '/QC/SummaryStats'
        }

        call mapping.mapping as indv_mapping {
            input :
                fastqfile=eachfastq,
                index_files=bowtie_index_,
                metricsfile=indv_bfs.metrics_out,
                blacklist=blacklist,
                default_location='SAMPLE/' + sub(basename(eachfastq),'_R?[12].*\.f.*q\.gz','') + '/BAM_files'
        }

        call fastqc.fastqc as indv_bamfqc {
            input :
                inputfile=indv_mapping.sorted_bam,
                default_location='SAMPLE/' + sub(basename(eachfastq),'_R?[12].*\.f.*q\.gz','') + '/QC/FastQC'
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
                default_location='SAMPLE/' + sub(basename(eachfastq),'_R?[12].*\.f.*q\.gz','') + '/QC/SummaryStats'
        }
    } # end scatter (for eachfastq)

    # MERGE BAM files
    # Execute analysis on merge bam file
    # Analysis executed:
    #   Merge BAM (SE mode : for each fastq paired-end)
    #   FastQC on Merge BAM (AllMerge_<number>_mapped)

    call util.mergehtml {
        input:
            htmlfiles=indv_summarystats.xhtml,
            txtfiles=indv_summarystats.textfile,
            default_location='SAMPLE',
            outputfile = 'AllMapped_peaseq-summary-stats.html'
    }

    call samtools.mergebam as SE_mergebam {
        input:
            bamfiles=indv_mapping.sorted_bam,
            default_location = if defined(results_name) then results_name + '_SE/BAM_files' else 'AllMapped_' + length(indv_mapping.sorted_bam) + '_SE/BAM_files',
            outputfile = if defined(results_name) then results_name + '_SE.sorted.bam' else 'AllMapped_' + length(all_sample_fastqfiles) + '_SE.sorted.bam'
    }

    call fastqc.fastqc as SE_mergebamfqc {
        input:
            inputfile=SE_mergebam.mergebam,
            default_location=sub(basename(SE_mergebam.mergebam),'\.sorted\.b.*$','') + '/QC/FastQC'
    }

    call samtools.indexstats as SE_mergeindexstats {
        input:
            bamfile=SE_mergebam.mergebam,
            default_location=sub(basename(SE_mergebam.mergebam),'\.sorted\.b.*$','') + '/BAM_files'
    }

    if ( defined(blacklist) ) {
        # remove blacklist regions
        String string_blacklist = "" #buffer to allow for blacklist optionality
        File blacklist_ = select_first([blacklist, string_blacklist])
        call bedtools.intersect as SE_merge_rmblklist {
            input :
                fileA=SE_mergebam.mergebam,
                fileB=blacklist_,
                default_location=sub(basename(SE_mergebam.mergebam),'\.sorted\.b.*$','') + '/BAM_files',
                        nooverlap=true
        }
        call samtools.indexstats as SE_merge_bklist {
            input :
                bamfile=SE_merge_rmblklist.intersect_out,
                default_location=sub(basename(SE_mergebam.mergebam),'\.sorted\.b.*$','') + '/BAM_files'
        }
    } # end if blacklist provided

    File SE_mergebam_afterbklist = select_first([SE_merge_rmblklist.intersect_out, SE_mergebam.mergebam])

    call samtools.markdup as SE_merge_markdup {
        input :
            bamfile=SE_mergebam_afterbklist,
            default_location=sub(basename(SE_mergebam_afterbklist),'\.sorted\.b.*$','') + '/BAM_files'
        }

    call samtools.indexstats as SE_merge_mkdup {
        input :
            bamfile=SE_merge_markdup.mkdupbam,
            default_location=sub(basename(SE_mergebam_afterbklist),'\.sorted\.b.*$','') + '/BAM_files'
    }

### ------------------------------------------------- ###
### ---------------- S E C T I O N 3 ---------------- ###
### ---- Paired End (PE) Mode for all fastqfiles ---- ###
### ---- A: analysis if multiple FASTQs provided ---- ###
### ------------------------------------------------- ###

    # if multiple fastqfiles are provided
    Boolean multi_fastqpair = if length(sample_fastqfiles) > 1 then true else false
    Boolean one_fastqpair = if length(sample_fastqfiles) == 1 then true else false

############# multi fastq

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

        call mapping.mapping as uno_PE_mapping {
            input :
                fastqfile=sample_fastqfiles[0].left,
                fastqfile_R2=sample_fastqfiles[0].right,
                insert_size=bowtie_insert_size,
                index_files=bowtie_index_,
                blacklist=blacklist,
                paired_end=true,
                results_name=results_name,
                default_location=if defined(results_name) then results_name + '_PE/BAM_files' else sub(basename(sample_fastqfiles[0].left),'_R?[12].*\.f.*q\.gz','') + '_PE/BAM_files'
        }

        call fastqc.fastqc as uno_PE_bamfqc {
            input :
                inputfile=uno_PE_mapping.sorted_bam,
                default_location=sub(basename(uno_PE_mapping.sorted_bam),'\.sorted\.bam','') + '/QC/FastQC'
        }

        call runspp.runspp as uno_PE_runspp {
            input:
                bamfile=select_first([uno_PE_mapping.bklist_bam, uno_PE_mapping.sorted_bam])
        }

        call bedtools.bamtobed as uno_PE_bamtobed {
            input:
                bamfile=select_first([uno_PE_mapping.bklist_bam, uno_PE_mapping.sorted_bam])
        }
    } # end if length(fastqfiles) == 1: one_fastqpair

### ------------------------------------------------- ###
### ---------------- S E C T I O N 4 ---------------- ###
### --------------- ChIP-seq Analysis --------------- ###
### ------------------------------------------------- ###



### ------------------------------------------------- ###
### ---------------- S E C T I O N 5 ---------------- ###
### ------------------ OUTPUT FILES ----------------- ###
### ------------------------------------------------- ###

    output {
        #FASTQC
        Array[File?]? indv_s_fastqc = indv_fastqc.htmlfile
        Array[File?]? indv_s_zipfile = indv_fastqc.zipfile
        Array[File?]? indv_s_bam_htmlfile = indv_bamfqc.htmlfile
        Array[File?]? indv_s_bam_zipfile = indv_bamfqc.zipfile
        File? s_mergebam_htmlfile = SE_mergebamfqc.htmlfile
        File? s_mergebam_zipfile = SE_mergebamfqc.zipfile

        File? uno_sp_bam_htmlfile = uno_PE_bamfqc.htmlfile
        File? uno_sp_bam_zipfile = uno_PE_bamfqc.zipfile

        #BASICMETRICS
        Array[File?]? s_metrics_out = indv_bfs.metrics_out

        #BAMFILES
        Array[File?]? indv_s_sortedbam = indv_mapping.sorted_bam
        Array[File?]? indv_s_indexbam = indv_mapping.bam_index
        Array[File?]? indv_s_bkbam = indv_mapping.bklist_bam
        Array[File?]? indv_s_bkindexbam = indv_mapping.bklist_index
        Array[File?]? indv_s_rmbam = indv_mapping.mkdup_bam
        Array[File?]? indv_s_rmindexbam = indv_mapping.mkdup_index

        File? uno_s_sortedbam = uno_PE_mapping.sorted_bam
        File? uno_s_indexbam = uno_PE_mapping.bam_index
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

        #QC-STATS
        Array[File?]? s_qc_statsfile = indv_summarystats.statsfile
        Array[File?]? s_qc_htmlfile = indv_summarystats.htmlfile
        Array[File?]? s_qc_textfile = indv_summarystats.textfile
        File? s_qc_mergehtml = mergehtml.mergefile

    }
}
