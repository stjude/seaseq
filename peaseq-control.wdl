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
        File? spikein_reference
        File? blacklist
        File gtf
        Array[File]? bowtie_index
        Array[File]? spikein_bowtie_index
        Array[File]? motif_databases

        # group: input_genomic_data
        Array[String]? sample_sraid
        Array[File]? sample_R1_fastq
        Array[File]? sample_R2_fastq
        Array[String]? control_sraid
        Array[File]? control_R1_fastq
        Array[File]? control_R2_fastq

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
        control_sraid: {
            description: 'One or more input/control SRA (Sequence Read Archive) run identifiers',
            group: 'input_genomic_data',
            help: 'Input publicly available FASTQs (SRRs). Multiple SRRs are separated by commas (,).',
            example: 'SRR12345678'
        }
        control_R1_fastq: {
            description: 'One or more input/control R1 FASTQs',
            group: 'input_genomic_data',
            help: 'Upload zipped FASTQ files.',
            patterns: ["*.fq.gz", "*.fastq.gz"]
        }
        control_R2_fastq: {
            description: 'One or more input/control R2 FASTQs',
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

    if ( defined(control_sraid) ) {
        # Download control file(s) from SRA database
        # outputs:
        #    fastqdump.fastqfile : downloaded sample files in fastq.gz format
        Array[String] c_sra = [1] #buffer to allow for sra_id optionality
        Array[String] c_sraid = select_first([control_sraid, c_sra])
        scatter (eachsra in c_sraid) {
            call sra.fastqdump as c_fastqdump {
                input :
                    sra_id=eachsra,
                    cloud=false
            }
            File c_R1end = select_first([c_fastqdump.R1end, c_sra[0]])
            File c_R2end = select_first([c_fastqdump.R2end, c_sra[0]])
        } # end scatter each sra

        Array[File] control_R1_srafile = c_R1end
        Array[File] control_R2_srafile = c_R2end
    } # end if control_sra_id

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

    # Spike-in DNA
    #3. Bowtie INDEX files if not provided
    String string_spikein = "1"
    Array[String] string_spikein_buffer = [1]
    if ( !defined(spikein_bowtie_index) && defined(spikein_reference) ) {
        # create bowtie index on spikein genome
        call bowtie.index as spikein_bowtie_idx {
            input :
                reference=select_first([spikein_reference, string_spikein])
        }
    }

    #4. Make sure indexes are six else build indexes for Spike-in DNA
    if ( defined(spikein_bowtie_index) ) {
        # check total number of bowtie indexes provided
        Array[File] int_spikein_bowtie_index = select_first([spikein_bowtie_index, string_spikein_buffer])
        if ( length(int_spikein_bowtie_index) != 6 ) {
            # create bowtie index if 6 index files aren't provided
            call bowtie.index as spikein_bowtie_idx_2 {
                input :
                    reference=select_first([spikein_reference, string_spikein])
            }
        }
    }
    Array[File] actual_spikein_bowtie_index = select_first([spikein_bowtie_idx_2.bowtie_indexes, spikein_bowtie_idx.bowtie_indexes, spikein_bowtie_index, string_spikein_buffer])

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

        # Order FASTQs
        if (length(sample_R1_fastqfile) > 1) {
            call peaseq_util.sortfiles as R1_sorted { input: fastqfiles=sample_R1_fastqfile }
            call peaseq_util.sortfiles as R2_sorted { input: fastqfiles=sample_R2_fastqfile }
        }
        Array[File] sample_R1_fastqfiles = select_first([R1_sorted.allfiles, sample_R1_fastqfile])
        Array[File] sample_R2_fastqfiles = select_first([R2_sorted.allfiles, sample_R2_fastqfile])
    }

    if ( defined(control_R1_fastq) ) {
        Array[String] c_string_fastq = [1] #buffer to allow for fastq optionality
        Array[File] c_R1_fastq = select_first([control_R1_fastq, c_string_fastq])
        Array[File] control_R1_fastqfile = c_R1_fastq
        Array[File] c_R2_fastq = select_first([control_R2_fastq, c_string_fastq])
        Array[File] control_R2_fastqfile = c_R2_fastq

        # Order FASTQs
        if (length(control_R1_fastqfile) > 1) {
            call peaseq_util.sortfiles as c_R1_sorted { input: fastqfiles=control_R1_fastqfile }
            call peaseq_util.sortfiles as c_R2_sorted { input: fastqfiles=control_R2_fastqfile }
        }
        Array[File] control_R1_fastqfiles = select_first([c_R1_sorted.allfiles, control_R1_fastqfile])
        Array[File] control_R2_fastqfiles = select_first([c_R2_sorted.allfiles, control_R2_fastqfile])
    }

    # collate all fastqfiles
    Array[File] original_sample_R1 = flatten(select_all([sample_R1_srafile, sample_R1_fastqfiles]))
    Array[File] original_sample_R2 = flatten(select_all([sample_R2_srafile, sample_R2_fastqfiles]))
    Array[File] original_all_sample_fastqfiles = flatten(select_all([sample_R1_srafile, sample_R1_fastqfiles,sample_R2_srafile, sample_R2_fastqfiles]))
    Array[File] original_control_R1 = flatten(select_all([control_R1_srafile, control_R1_fastqfiles]))
    Array[File] original_control_R2 = flatten(select_all([control_R2_srafile, control_R2_fastqfiles]))
    Array[File] original_all_control_fastqfiles = flatten(select_all([control_R1_srafile, control_R1_fastqfiles,control_R2_srafile, control_R2_fastqfiles]))

    # transpose to paired-end tuples
    Array[Pair[File, File]] original_sample_fastqfiles = zip(original_sample_R1, original_sample_R2)
    Array[Pair[File, File]] original_control_fastqfiles = zip(original_control_R1, original_control_R2)

### ------------------------------------------------- ###
### ---------------- S E C T I O N 1 ---------------- ###
### ----------- B: remove Spike-IN reads ------------ ###
### ------------------------------------------------- ###

    # if multiple fastqfiles are provided
    Boolean multi_fastqpair = if length(original_sample_fastqfiles) > 1 then true else false
    Boolean one_fastqpair = if length(original_sample_fastqfiles) == 1 then true else false
    Boolean multi_control_fastqpair = if length(original_control_fastqfiles) > 1 then true else false
    Boolean one_control_fastqpair = if length(original_control_fastqfiles) == 1 then true else false

    if ( defined(spikein_bowtie_index) || defined(spikein_reference) ) {
        scatter (eachfastq in original_all_sample_fastqfiles) {
            call fastqc.fastqc as spikein_s_indv_fastqc {
                input :
                    inputfile=eachfastq,
                    default_location=if multi_fastqpair then 'SAMPLE/individual_fastqs/' + sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/SpikeIn/QC/FastQC' else 'SAMPLE/' + sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/SpikeIn/QC/FastQC'
            }
        }
        scatter (eachfastq in original_all_control_fastqfiles) {
            call fastqc.fastqc as spikein_c_indv_fastqc {
                input :
                    inputfile=eachfastq,
                    default_location=if multi_control_fastqpair then 'CONTROL/individual_fastqs/' + sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/SpikeIn/QC/FastQC' else 'CONTROL/' + sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/SpikeIn/QC/FastQC'
            }
        }
        scatter (fastqpair in original_sample_fastqfiles) {
            call util.basicfastqstats as spikein_s_indv_R1_bfs {
                input :
                    fastqfile=fastqpair.left,
                    default_location=if multi_fastqpair then 'SAMPLE/individual_fastqs/' + sub(basename(fastqpair.left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/SpikeIn/QC/SummaryStats' else 'SAMPLE/' + sub(basename(fastqpair.left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/SpikeIn/QC/SummaryStats'
            }
            call bowtie.spikein_PE as spikein_s_indv_map {
                input :
                    fastqfile=fastqpair.left,
                    fastqfile_R2=fastqpair.right,
                    metricsfile=spikein_s_indv_R1_bfs.metrics_out,
                    insert_size=insertsize,
                    strandedness=strandedness,
                    index_files=actual_spikein_bowtie_index,
                    default_location=if multi_fastqpair then 'SAMPLE/individual_fastqs/' + sub(basename(fastqpair.left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/SpikeIn/QC/SummaryStats' else 'SAMPLE/' + sub(basename(fastqpair.left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/SpikeIn/QC/SummaryStats'
            }
        }
        scatter (fastqpair in original_control_fastqfiles) {
            call util.basicfastqstats as spikein_c_indv_R1_bfs {
                input :
                    fastqfile=fastqpair.left,
                    default_location=if multi_control_fastqpair then 'CONTROL/individual_fastqs/' + sub(basename(fastqpair.left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/SpikeIn/QC/SummaryStats' else 'CONTROL/' + sub(basename(fastqpair.left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/SpikeIn/QC/SummaryStats'
            }
            call bowtie.spikein_PE as spikein_c_indv_map {
                input :
                    fastqfile=fastqpair.left,
                    fastqfile_R2=fastqpair.right,
                    metricsfile=spikein_c_indv_R1_bfs.metrics_out,
                    insert_size=insertsize,
                    strandedness=strandedness,
                    index_files=actual_spikein_bowtie_index,
                    default_location=if multi_control_fastqpair then 'CONTROL/individual_fastqs/' + sub(basename(fastqpair.left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/SpikeIn/QC/SummaryStats' else 'CONTROL/' + sub(basename(fastqpair.left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/SpikeIn/QC/SummaryStats'
            }
        }

        if (length(original_sample_R1) > 1) {
            call peaseq_util.sortfiles as spikein_s_R1_sorted { input: fastqfiles=spikein_s_indv_map.unaligned_1 }
            call peaseq_util.sortfiles as spikein_s_R2_sorted { input: fastqfiles=spikein_s_indv_map.unaligned_2 }
        }
        if (length(original_control_R1) > 1) {
            call peaseq_util.sortfiles as spikein_c_R1_sorted { input: fastqfiles=spikein_c_indv_map.unaligned_1 }
            call peaseq_util.sortfiles as spikein_c_R2_sorted { input: fastqfiles=spikein_c_indv_map.unaligned_2 }
        }

        Array[File] spikein_sample_R1 = select_first([spikein_s_R1_sorted.allfiles, spikein_s_indv_map.unaligned_1])
        Array[File] spikein_sample_R2 = select_first([spikein_s_R2_sorted.allfiles, spikein_s_indv_map.unaligned_2])
        Array[File] spikein_control_R1 = select_first([spikein_c_R1_sorted.allfiles, spikein_c_indv_map.unaligned_1])
        Array[File] spikein_control_R2 = select_first([spikein_c_R2_sorted.allfiles, spikein_c_indv_map.unaligned_2])

        Array[File] spikein_all_sample_fastqfiles = flatten(select_all([spikein_sample_R1, spikein_sample_R2]))
        Array[Pair[File, File]] spikein_sample_fastqfiles = zip(spikein_sample_R1, spikein_sample_R2)
        Array[File] spikein_all_control_fastqfiles = flatten(select_all([spikein_control_R1, spikein_control_R2]))
        Array[Pair[File, File]] spikein_control_fastqfiles = zip(spikein_control_R1, spikein_control_R2)
    }
        
    Array[File] all_sample_fastqfiles = select_first([spikein_all_sample_fastqfiles, original_all_sample_fastqfiles])
    Array[Pair[File, File]] sample_fastqfiles = select_first([spikein_sample_fastqfiles, original_sample_fastqfiles])
    Array[File] all_control_fastqfiles = select_first([spikein_all_control_fastqfiles, original_all_control_fastqfiles])
    Array[Pair[File, File]] control_fastqfiles = select_first([spikein_control_fastqfiles, original_control_fastqfiles])

### ------------------------------------------------- ###
### ---------------- S E C T I O N 2 ---------------- ###
### ---- Single End (SE) Mode for all fastqfiles ---- ###
### ------------------------------------------------- ###

    scatter (eachfastq in all_sample_fastqfiles) {
        call fastqc.fastqc as s_indv_fastqc {
            input :
                inputfile=eachfastq,
                default_location=if multi_fastqpair then 'SAMPLE/individual_fastqs/' + sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/FastQC' else 'SAMPLE/' + sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/FastQC'
        }

        call util.basicfastqstats as indv_s_bfs {
            input :
                fastqfile=eachfastq,
                default_location=if multi_fastqpair then 'SAMPLE/individual_fastqs/' + sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/SummaryStats' else 'SAMPLE/' + sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/SummaryStats'
        }

        call mapping.mapping as indv_s_mapping {
            input :
                fastqfile=eachfastq,
                index_files=actual_bowtie_index,
                metricsfile=indv_s_bfs.metrics_out,
                blacklist=blacklist,
                default_location=if multi_fastqpair then  'SAMPLE/individual_fastqs/' + sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/BAM_files' else 'SAMPLE/' + sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/BAM_files'
        }

        call fastqc.fastqc as indv_s_bamfqc {
            input :
                inputfile=indv_s_mapping.sorted_bam,
                default_location=if multi_fastqpair then  'SAMPLE/individual_fastqs/' + sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/FastQC' else 'SAMPLE/' + sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/FastQC'
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
                fastq_type="PEAseq Sample FASTQ",
                bambed=indv_s_bamtobed.bedfile,
                sppfile=indv_s_runspp.spp_out,
                fastqczip=s_indv_fastqc.zipfile,
                bamflag=indv_s_mapping.bam_stats,
                rmdupflag=indv_s_mapping.mkdup_stats,
                bkflag=indv_s_mapping.bklist_stats,
                fastqmetrics=indv_s_bfs.metrics_out,
                default_location=if multi_fastqpair then 'SAMPLE/individual_fastqs/' + sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/QC/SummaryStats' else 'SAMPLE/' + sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/QC/SummaryStats'
        }
    } # end scatter (for eachfastq : samples)

    # MERGE BAM files
    # Execute analysis on merge bam file
    # Analysis executed:
    #   Merge BAM (SE mode : for each fastq paired-end)

    call util.mergehtml as SE_s_mergehtml {
        input:
            htmlfiles=indv_s_summarystats.xhtml,
            txtfiles=indv_s_summarystats.textfile,
            default_location='SAMPLE',
            outputfile = 'AllCases-summary-stats.html'
    }

    call samtools.mergebam as SE_s_mergebam {
        input:
            bamfiles=indv_s_mapping.sorted_bam,
            metricsfiles=indv_s_bfs.metrics_out,
            default_location=if multi_fastqpair then 'SAMPLE/' + 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs/single-end_mode/BAM_files' else 'SAMPLE/' + sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/BAM_files',
            outputfile='AllCases_' + length(all_sample_fastqfiles) + 'fastqs_SE.sorted.bam'
    }

    call fastqc.fastqc as SE_s_mergebamfqc {
        input:
            inputfile=SE_s_mergebam.mergebam,
            default_location= if multi_fastqpair then 'SAMPLE/' + 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs/single-end_mode/QC/FastQC' else 'SAMPLE/' + sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/QC/FastQC'
    }

    call samtools.indexstats as SE_s_mergeindexstats {
        input:
            bamfile=SE_s_mergebam.mergebam,
            default_location=if multi_fastqpair then 'SAMPLE/' + 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs/single-end_mode/BAM_files' else 'SAMPLE/' + sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/BAM_files'
    }

    if ( defined(blacklist) ) {
        # remove blacklist regions
        String string_blacklist = "" #buffer to allow for blacklist optionality
        File blacklist_file = select_first([blacklist, string_blacklist])
        call bedtools.intersect as SE_s_merge_rmblklist {
            input :
                fileA=SE_s_mergebam.mergebam,
                fileB=blacklist_file,
                default_location=if multi_fastqpair then 'SAMPLE/' + 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs/single-end_mode/BAM_files' else 'SAMPLE/' + sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/BAM_files',
                nooverlap=true
        }
        call samtools.indexstats as SE_s_merge_bklist {
            input :
                bamfile=SE_s_merge_rmblklist.intersect_out,
                default_location=if multi_fastqpair then 'SAMPLE/' + 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs/single-end_mode/BAM_files' else 'SAMPLE/' + sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/BAM_files'
        }
    } # end if blacklist provided

    File SE_s_mergebam_afterbklist = select_first([SE_s_merge_rmblklist.intersect_out, SE_s_mergebam.mergebam])

    call samtools.markdup as SE_s_merge_markdup {
        input :
            bamfile=SE_s_mergebam_afterbklist,
            default_location=if multi_fastqpair then 'SAMPLE/' + 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs/single-end_mode/BAM_files' else 'SAMPLE/' + sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/BAM_files'
    }

    call samtools.indexstats as SE_s_merge_mkdup {
        input :
            bamfile=SE_s_merge_markdup.mkdupbam,
            default_location=if multi_fastqpair then 'SAMPLE/' + 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs/single-end_mode/BAM_files' else 'SAMPLE/' + sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/BAM_files'
    }

    # CONTROL FASTQ FILES
    scatter (eachfastq in all_control_fastqfiles) {
        call fastqc.fastqc as c_indv_fastqc {
            input :
                inputfile=eachfastq,
                default_location=if multi_control_fastqpair then 'CONTROL/individual_fastqs/' + sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/FastQC' else 'CONTROL/' + sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/FastQC'
        }

        call util.basicfastqstats as indv_c_bfs {
            input :
                fastqfile=eachfastq,
                default_location=if multi_control_fastqpair then 'CONTROL/individual_fastqs/' + sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/SummaryStats' else 'CONTROL/' + sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/SummaryStats'
        }

        call mapping.mapping as indv_c_mapping {
            input :
                fastqfile=eachfastq,
                index_files=actual_bowtie_index,
                metricsfile=indv_c_bfs.metrics_out,
                blacklist=blacklist,
                default_location=if multi_control_fastqpair then 'CONTROL/individual_fastqs/' + sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/BAM_files' else 'CONTROL/' + sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/BAM_files'
        }

        call fastqc.fastqc as indv_c_bamfqc {
            input :
                inputfile=indv_c_mapping.sorted_bam,
                default_location=if multi_control_fastqpair then 'CONTROL/individual_fastqs/' + sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/FastQC' else 'CONTROL/' + sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/FastQC'
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
                fastq_type="PEAseq Control FASTQ",
                bambed=indv_c_bamtobed.bedfile,
                sppfile=indv_c_runspp.spp_out,
                fastqczip=c_indv_fastqc.zipfile,
                bamflag=indv_c_mapping.bam_stats,
                rmdupflag=indv_c_mapping.mkdup_stats,
                bkflag=indv_c_mapping.bklist_stats,
                fastqmetrics=indv_c_bfs.metrics_out,
                default_location=if multi_control_fastqpair then 'CONTROL/individual_fastqs/' + sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/QC/SummaryStats' else 'CONTROL/' + sub(basename(eachfastq),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/QC/SummaryStats'
        }
    } # end scatter (for eachfastq : control)

    # MERGE BAM files
    # Execute analysis on merge bam file
    # Analysis executed:
    #   Merge BAM (SE mode : for each fastq paired-end)

    call util.mergehtml as SE_c_mergehtml {
        input:
            fastq_type="PEAseq Control FASTQs",
            htmlfiles=indv_c_summarystats.xhtml,
            txtfiles=indv_c_summarystats.textfile,
            default_location='CONTROL',
            outputfile = 'AllControls-summary-stats.html'
    }

    call samtools.mergebam as SE_c_mergebam {
        input:
            bamfiles=indv_c_mapping.sorted_bam,
            metricsfiles=indv_c_bfs.metrics_out,
            default_location=if multi_control_fastqpair then 'CONTROL/' + 'AllControls_' + length(control_fastqfiles) + 'fastqpairs/single-end_mode/BAM_files' else 'CONTROL/' + sub(basename(control_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/BAM_files',
            outputfile='AllControls_' + length(all_control_fastqfiles) + 'fastqs_SE.sorted.bam'
    }

    call fastqc.fastqc as SE_c_mergebamfqc {
        input:
            inputfile=SE_c_mergebam.mergebam,
            default_location=if multi_control_fastqpair then 'CONTROL/' + 'AllControls_' + length(control_fastqfiles) + 'fastqpairs/single-end_mode/QC/FastQC' else 'CONTROL/' + sub(basename(control_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/QC/FastQC'
    }

    call samtools.indexstats as SE_c_mergeindexstats {
        input:
            bamfile=SE_c_mergebam.mergebam,
            default_location=if multi_control_fastqpair then 'CONTROL/' + 'AllControls_' + length(control_fastqfiles) + 'fastqpairs/single-end_mode/BAM_files' else 'CONTROL/' + sub(basename(control_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/BAM_files'
    }

    if ( defined(blacklist) ) {
        # remove blacklist regions
        String c_string_blacklist = "" #buffer to allow for blacklist optionality
        File c_blacklist_file = select_first([blacklist, c_string_blacklist])
        call bedtools.intersect as SE_c_merge_rmblklist {
            input :
                fileA=SE_c_mergebam.mergebam,
                fileB=c_blacklist_file,
                default_location=if multi_control_fastqpair then 'CONTROL/' + 'AllControls_' + length(control_fastqfiles) + 'fastqpairs/single-end_mode/BAM_files' else 'CONTROL/' + sub(basename(control_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/BAM_files',
                nooverlap=true
        }
        call samtools.indexstats as SE_c_merge_bklist {
            input :
                bamfile=SE_c_merge_rmblklist.intersect_out,
                default_location=if multi_control_fastqpair then 'CONTROL/' + 'AllControls_' + length(control_fastqfiles) + 'fastqpairs/single-end_mode/BAM_files' else 'CONTROL/' + sub(basename(control_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/BAM_files'
        }
    } # end if blacklist provided

    File SE_c_mergebam_afterbklist = select_first([SE_c_merge_rmblklist.intersect_out, SE_c_mergebam.mergebam])

    call samtools.markdup as SE_c_merge_markdup {
        input :
            bamfile=SE_c_mergebam_afterbklist,
            default_location=if multi_control_fastqpair then 'CONTROL/' + 'AllControls_' + length(control_fastqfiles) + 'fastqpairs/single-end_mode/BAM_files' else 'CONTROL/' + sub(basename(control_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/BAM_files'
        }

    call samtools.indexstats as SE_c_merge_mkdup {
        input :
            bamfile=SE_c_merge_markdup.mkdupbam,
            default_location=if multi_control_fastqpair then 'CONTROL/' + 'AllControls_' + length(control_fastqfiles) + 'fastqpairs/single-end_mode/BAM_files' else 'CONTROL/' + sub(basename(control_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/BAM_files'
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

            call util.basicfastqstats as s_indv_R1_bfs {
                input :
                    fastqfile=fastqpair.left,
                    default_location=if multi_fastqpair then 'SAMPLE/individual_fastqs/' + sub(basename(fastqpair.left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/SummaryStats' else 'SAMPLE/' + sub(basename(fastqpair.left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/SummaryStats'
            }

            call mapping.mapping as s_indv_PE_mapping {
                input :
                    fastqfile=fastqpair.left,
                    fastqfile_R2=fastqpair.right,
                    metricsfile=s_indv_R1_bfs.metrics_out,
                    insert_size=insertsize,
                    strandedness=strandedness,
                    index_files=actual_bowtie_index,
                    blacklist=blacklist,
                    paired_end=true,
                    default_location=if multi_fastqpair then 'SAMPLE/individual_fastqs/' + sub(basename(fastqpair.left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/BAM_files' else 'SAMPLE/' + sub(basename(fastqpair.left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/BAM_files'
            }

            call fastqc.fastqc as s_indv_PE_bamfqc {
                input :
                    inputfile=s_indv_PE_mapping.sorted_bam,
                    default_location=if multi_fastqpair then 'SAMPLE/individual_fastqs/' + sub(basename(fastqpair.left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/FastQC' else 'SAMPLE/' + sub(basename(fastqpair.left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/FastQC'
            }

            call runspp.runspp as s_indv_PE_runspp {
                input:
                    bamfile=select_first([s_indv_PE_mapping.bklist_bam, s_indv_PE_mapping.sorted_bam])
            }

            call peaseq_util.pe_bamtobed as s_indv_PE_bamtobed {
                input:
                    bamfile=select_first([s_indv_PE_mapping.bklist_bam, s_indv_PE_mapping.sorted_bam])
            }

            call util.evalstats as s_indv_PE_summarystats {
                input:
                    fastq_type="PEAseq Sample FASTQ",
                    bambed=s_indv_PE_bamtobed.bedfile,
                    sppfile=s_indv_PE_runspp.spp_out,
                    fastqczip=s_indv_PE_bamfqc.zipfile,
                    bamflag=s_indv_PE_mapping.bam_stats,
                    rmdupflag=s_indv_PE_mapping.mkdup_stats,
                    bkflag=s_indv_PE_mapping.bklist_stats,
                    default_location=if multi_fastqpair then 'SAMPLE/individual_fastqs/' + sub(basename(fastqpair.left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/SummaryStats' else 'SAMPLE/' + sub(basename(fastqpair.left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/SummaryStats'
            }
        } # end scatter (for each sample fastq)

        # MERGE BAM FILES
        # Execute analysis on merge bam file
        # Analysis executed:
        #   Merge BAM (if more than 1 fastq is provided)
        #   FastQC on Merge BAM (AllCases_<number>_mapped)

        # merge bam files and perform fasTQC if more than one is provided
        call util.mergehtml as s_PE_mergehtml {
            input:
                htmlfiles=s_indv_PE_summarystats.xhtml,
                txtfiles=s_indv_PE_summarystats.textfile,
                default_location='SAMPLE',
                outputfile = 'AllCases_PEmode_peaseq-summary-stats.html'
        }

        call peaseq_util.pe_mergehtml as s_final_mergehtml {
            input:
                pe_htmlfiles=s_indv_PE_summarystats.xhtml,
                pe_txtfiles=s_indv_PE_summarystats.textfile,
                se_htmlfiles=indv_s_summarystats.xhtml,
                se_txtfiles=indv_s_summarystats.textfile,
                default_location='SAMPLE',
                outputfile = 'AllCases-summary-stats.html'
        }

        call samtools.mergebam as s_PE_mergebam {
            input:
                bamfiles=s_indv_PE_mapping.as_sortedbam,
                metricsfiles=indv_s_bfs.metrics_out,
                paired_end=true,
                default_location='SAMPLE/' + 'AllCases_' + length(s_indv_PE_mapping.sorted_bam) + 'fastqpairs/BAM_files',
                outputfile='AllCases_' + length(sample_fastqfiles) + 'fastqpairs.sorted.bam',
                fixmatefile='AllCases_' + length(sample_fastqfiles) + 'fastqpairs.fixmate.bam'
        }

        call fastqc.fastqc as s_PE_mergebamfqc {
            input:
	        inputfile=s_PE_mergebam.mergebam,
                default_location='SAMPLE/' + sub(basename(s_PE_mergebam.mergebam),'.sorted.b.*$','') + '/QC/FastQC'
        }

        call samtools.indexstats as s_PE_mergeindexstats {
            input:
                bamfile=s_PE_mergebam.mergebam,
                default_location='SAMPLE/' + sub(basename(s_PE_mergebam.mergebam),'.sorted.b.*$','') + '/BAM_files'
        }

        if ( defined(blacklist) ) {
            # remove blacklist regions
            String s_string_pe_blacklist = "" #buffer to allow for blacklist optionality
            File s_blacklist_pe = select_first([blacklist, s_string_pe_blacklist])
            String s_string_pe_fixmate = "" #buffer to allow for blacklist optionality
            File s_fixmate_pe = select_first([s_PE_mergebam.fixmatemergebam, s_string_pe_fixmate])
            call bedtools.pairtobed as s_PE_merge_rmblklist {
                input :
                    fileA=s_fixmate_pe,
                    fileB=s_blacklist_pe,
                    default_location='SAMPLE/' + sub(basename(s_PE_mergebam.mergebam),'.sorted.b.*$','') + '/BAM_files'
            }
            call samtools.indexstats as s_PE_merge_bklist {
                input :
                    bamfile=s_PE_merge_rmblklist.pairtobed_out,
                    default_location='SAMPLE/' + sub(basename(s_PE_mergebam.mergebam),'.sorted.b.*$','') + '/BAM_files'
            }
        } # end if blacklist provided

        File s_PE_mergebam_afterbklist = select_first([s_PE_merge_rmblklist.pairtobed_out, s_PE_mergebam.mergebam])

        call samtools.markdup as s_PE_merge_markdup {
            input :
                bamfile=s_PE_mergebam_afterbklist,
                default_location='SAMPLE/' + sub(basename(s_PE_mergebam.mergebam),'.sorted.b.*$','') + '/BAM_files'
        }

        call samtools.indexstats as s_PE_merge_mkdup {
            input :
                bamfile=s_PE_merge_markdup.mkdupbam,
                default_location='SAMPLE/' + sub(basename(s_PE_mergebam.mergebam),'.sorted.b.*$','') + '/BAM_files'
        }
    } # end if length(fastqfiles) > 1: multi_fastqpair

    # CONTROL FASTQ FILES
    if ( multi_control_fastqpair ) {
        scatter (fastqpair in control_fastqfiles) {
            # Execute analysis on each fastq file provided
            # Analysis executed:
            #   Reference Alignment using Bowtie (-k2 -m2)
            #   Convert SAM to BAM
            #   FastQC on BAM files
            #   Remove Blacklists (if provided)
            #   Remove read duplicates
            #   Summary statistics on FASTQs
            #   Combine html files into one for easy viewing

            call util.basicfastqstats as c_indv_R1_bfs {
                input :
                    fastqfile=fastqpair.left,
                    default_location=if multi_fastqpair then 'CONTROL/individual_fastqs/' + sub(basename(fastqpair.left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/SummaryStats' else 'CONTROL/' + sub(basename(fastqpair.left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/SummaryStats'
            }
            
            call mapping.mapping as c_indv_PE_mapping {
                input :
                    fastqfile=fastqpair.left,
                    fastqfile_R2=fastqpair.right,
                    metricsfile=c_indv_R1_bfs.metrics_out,
                    insert_size=insertsize,
                    strandedness=strandedness,
                    index_files=actual_bowtie_index,
                    blacklist=blacklist,
                    paired_end=true,
                    default_location=if multi_control_fastqpair then 'CONTROL/individual_fastqs/' + sub(basename(fastqpair.left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/BAM_files' else 'CONTROL/' + sub(basename(fastqpair.left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/BAM_files'
            }

            call fastqc.fastqc as c_indv_PE_bamfqc {
                input :
                    inputfile=c_indv_PE_mapping.sorted_bam,
                    default_location=if multi_control_fastqpair then 'CONTROL/individual_fastqs/' + sub(basename(fastqpair.left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/FastQC' else 'CONTROL/' + sub(basename(fastqpair.left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/FastQC'
            }

            call runspp.runspp as c_indv_PE_runspp {
                input:
                    bamfile=select_first([c_indv_PE_mapping.bklist_bam, c_indv_PE_mapping.sorted_bam])
            }

            call peaseq_util.pe_bamtobed as c_indv_PE_bamtobed {
                input:
                    bamfile=select_first([c_indv_PE_mapping.bklist_bam, c_indv_PE_mapping.sorted_bam])
            }

            call util.evalstats as c_indv_PE_summarystats {
                input:
                    fastq_type="PEAseq Control FASTQ",
                    bambed=c_indv_PE_bamtobed.bedfile,
                    sppfile=c_indv_PE_runspp.spp_out,
                    fastqczip=c_indv_PE_bamfqc.zipfile,
                    bamflag=c_indv_PE_mapping.bam_stats,
                    rmdupflag=c_indv_PE_mapping.mkdup_stats,
                    bkflag=c_indv_PE_mapping.bklist_stats,
                    default_location=if multi_control_fastqpair then 'CONTROL/individual_fastqs/' + sub(basename(fastqpair.left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/SummaryStats' else 'CONTROL/' + sub(basename(fastqpair.left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/SummaryStats'
            }
        } # end scatter (for each sample fastq)

        # MERGE BAM FILES
        # Execute analysis on merge bam file
        # Analysis executed:
        #   Merge BAM (if more than 1 fastq is provided)
        #   FastQC on Merge BAM (AllControl_<number>_mapped)

        # merge bam files and perform fasTQC if more than one is provided
        call util.mergehtml as c_PE_mergehtml {
            input:
                htmlfiles=c_indv_PE_summarystats.xhtml,
                txtfiles=c_indv_PE_summarystats.textfile,
                default_location='CONTROL',
                outputfile = 'AllControl_PEmode_peaseq-summary-stats.html'
        }

        call peaseq_util.pe_mergehtml as c_final_mergehtml {
            input:
                fastq_type = "PEAseq Control FASTQs",
                pe_htmlfiles=c_indv_PE_summarystats.xhtml,
                pe_txtfiles=c_indv_PE_summarystats.textfile,
                se_htmlfiles=indv_c_summarystats.xhtml,
                se_txtfiles=indv_c_summarystats.textfile,
                default_location='CONTROL',
                outputfile = 'AllControls-summary-stats.html'
        }

        call samtools.mergebam as c_PE_mergebam {
            input:
                bamfiles=c_indv_PE_mapping.as_sortedbam,
                metricsfiles=indv_c_bfs.metrics_out,
                paired_end=true,
                default_location = 'CONTROL/' + 'AllControls_' + length(c_indv_PE_mapping.sorted_bam) + 'fastqpairs/BAM_files',
                outputfile = 'AllControls_' + length(control_fastqfiles) + 'fastqpairs.sorted.bam',
                fixmatefile = 'AllControls_' + length(control_fastqfiles) + 'fastqpairs.fixmate.bam'
        }

        call fastqc.fastqc as c_PE_mergebamfqc {
            input:
	            inputfile=c_PE_mergebam.mergebam,
                default_location='CONTROL/' + sub(basename(c_PE_mergebam.mergebam),'.sorted.b.*$','') + '/QC/FastQC'
        }

        call samtools.indexstats as c_PE_mergeindexstats {
            input:
                bamfile=c_PE_mergebam.mergebam,
                default_location='CONTROL/' + sub(basename(c_PE_mergebam.mergebam),'.sorted.b.*$','') + '/BAM_files'
        }

        if ( defined(blacklist) ) {
            # remove blacklist regions
            String c_string_pe_blacklist = "" #buffer to allow for blacklist optionality
            File c_blacklist_pe = select_first([blacklist, c_string_pe_blacklist])
            String c_string_pe_fixmate = "" #buffer to allow for blacklist optionality
            File c_fixmate_pe = select_first([c_PE_mergebam.fixmatemergebam, c_string_pe_fixmate])
            call bedtools.pairtobed as c_PE_merge_rmblklist {
                input :
                    fileA=c_fixmate_pe,
                    fileB=c_blacklist_pe,
                    default_location='CONTROL/' + sub(basename(c_PE_mergebam.mergebam),'.sorted.b.*$','') + '/BAM_files'
            }
            call samtools.indexstats as c_PE_merge_bklist {
                input :
                    bamfile=c_PE_merge_rmblklist.pairtobed_out,
                    default_location='CONTROL/' + sub(basename(c_PE_mergebam.mergebam),'.sorted.b.*$','') + '/BAM_files'
            }
        } # end if blacklist provided

        File c_PE_mergebam_afterbklist = select_first([c_PE_merge_rmblklist.pairtobed_out, c_PE_mergebam.mergebam])

        call samtools.markdup as c_PE_merge_markdup {
            input :
                bamfile=c_PE_mergebam_afterbklist,
                default_location='CONTROL/' + sub(basename(c_PE_mergebam.mergebam),'.sorted.b.*$','') + '/BAM_files'
        }

        call samtools.indexstats as c_PE_merge_mkdup {
            input :
                bamfile=c_PE_merge_markdup.mkdupbam,
                default_location='CONTROL/' + sub(basename(c_PE_mergebam.mergebam),'.sorted.b.*$','') + '/BAM_files'
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

        call util.basicfastqstats as s_R1_bfs {
            input :
                fastqfile=sample_fastqfiles[0].left,
                default_location=if multi_fastqpair then 'SAMPLE/individual_fastqs/' + sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/SummaryStats' else 'SAMPLE/' + sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/SummaryStats'
                    
        }

        call mapping.mapping as s_uno_PE_mapping {
            input :
                fastqfile=sample_fastqfiles[0].left,
                fastqfile_R2=sample_fastqfiles[0].right,
                metricsfile=s_R1_bfs.metrics_out,
                insert_size=insertsize,
                strandedness=strandedness,
                index_files=actual_bowtie_index,
                blacklist=blacklist,
                paired_end=true,
                default_location=if multi_fastqpair then 'SAMPLE/individual_fastqs/' + sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/BAM_files' else 'SAMPLE/' + sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/BAM_files'
        }

        call fastqc.fastqc as s_uno_PE_bamfqc {
            input :
                inputfile=s_uno_PE_mapping.sorted_bam,
                default_location=if multi_fastqpair then 'SAMPLE/individual_fastqs/' + sub(basename(s_uno_PE_mapping.sorted_bam),'.sorted.bam','') + '/QC/FastQC' else 'SAMPLE/' + sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/FastQC'
        }

        call peaseq_util.pe_bamtobed as s_uno_PE_bamtobed {
            input:
                bamfile=select_first([s_uno_PE_mapping.bklist_bam, s_uno_PE_mapping.sorted_bam])
        }
    } # end if length(fastqfiles) == 1: one_fastqpair

    # CONTROL FASTQ FILE
    if ( one_control_fastqpair ) {
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

        call util.basicfastqstats as c_R1_bfs {
            input :
                fastqfile=control_fastqfiles[0].left,
                default_location=if multi_control_fastqpair then 'CONTROL/individual_fastqs/' + sub(basename(control_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/SummaryStats' else 'CONTROL/' + sub(basename(control_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/SummaryStats'
        }
        
        call mapping.mapping as c_uno_PE_mapping {
            input :
                fastqfile=control_fastqfiles[0].left,
                fastqfile_R2=control_fastqfiles[0].right,
                metricsfile=c_R1_bfs.metrics_out,
                insert_size=insertsize,
                strandedness=strandedness,
                index_files=actual_bowtie_index,
                blacklist=blacklist,
                paired_end=true,
                default_location=if multi_control_fastqpair then 'CONTROL/individual_fastqs/' + sub(basename(control_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/BAM_files' else 'CONTROL/' + sub(basename(control_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/BAM_files'
        }

        call fastqc.fastqc as c_uno_PE_bamfqc {
            input :
                inputfile=c_uno_PE_mapping.sorted_bam,
                default_location=if multi_control_fastqpair then 'CONTROL/individual_fastqs/' + sub(basename(c_uno_PE_mapping.sorted_bam),'.sorted.bam','') + '/QC/FastQC' else 'CONTROL/' + sub(basename(c_uno_PE_mapping.sorted_bam),'.sorted.bam','') + '/QC/FastQC'
        }

        call peaseq_util.pe_bamtobed as c_uno_PE_bamtobed {
            input:
                bamfile=select_first([c_uno_PE_mapping.bklist_bam, c_uno_PE_mapping.sorted_bam])
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
            bamfile=SE_s_mergebam_afterbklist,
            control=SE_c_mergebam_afterbklist,
            pvalue="1e-9",
            keep_dup="auto",
            egs=egs.genomesize,
            output_name=basename(SE_s_mergebam_afterbklist,'.bam') + '+control-p9_kd-auto',
            default_location = if defined(results_name) then results_name + '/single-end_mode/PEAKS/NARROW_peaks' + '/' + basename(SE_s_mergebam_afterbklist,'.bam') + '+control-p9_kd-auto' else if multi_fastqpair then 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs+control/single-end_mode/PEAKS/NARROW_peaks' + '/' + basename(SE_s_mergebam_afterbklist,'.bam') + '+control-p9_kd-auto' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '+control/single-end_mode/PEAKS/NARROW_peaks' + '/' + basename(SE_s_mergebam_afterbklist,'.bam') + '+control-p9_kd-auto',
            coverage_location = if defined(results_name) then results_name + '/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + basename(SE_s_mergebam_afterbklist,'.bam') + '+control-p9_kd-auto' else if multi_fastqpair then 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs+control/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + basename(SE_s_mergebam_afterbklist,'.bam') + '+control-p9_kd-auto' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '+control/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + basename(SE_s_mergebam_afterbklist,'.bam') + '+control-p9_kd-auto'
    }

    call util.addreadme as SE_addreadme {
        input :
            default_location = if defined(results_name) then results_name + '/single-end_mode/PEAKS' else if multi_fastqpair then 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs+control/single-end_mode/PEAKS' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '+control/single-end_mode/PEAKS'
    }

    call macs.macs as SE_all {
        input :
            bamfile=SE_s_mergebam_afterbklist,
            control=SE_c_mergebam_afterbklist,
            pvalue="1e-9",
            keep_dup="all",
            egs=egs.genomesize,
            output_name=basename(SE_s_mergebam_afterbklist,'.bam') + '+control-p9_kd-all',
            default_location = if defined(results_name) then results_name + '/single-end_mode/PEAKS/NARROW_peaks' + '/' + basename(SE_s_mergebam_afterbklist,'.bam') + '+control-p9_kd-all' else if multi_fastqpair then 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs+control/single-end_mode/PEAKS/NARROW_peaks' + '/' + basename(SE_s_mergebam_afterbklist,'.bam') + '+control-p9_kd-all' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '+control/single-end_mode/PEAKS/NARROW_peaks' + '/' + basename(SE_s_mergebam_afterbklist,'.bam') + '+control-p9_kd-all',
            coverage_location = if defined(results_name) then results_name + '/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + basename(SE_s_mergebam_afterbklist,'.bam') + '+control-p9_kd-all' else if multi_fastqpair then 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs+control/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + basename(SE_s_mergebam_afterbklist,'.bam') + '+control-p9_kd-all' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '+control/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + basename(SE_s_mergebam_afterbklist,'.bam') + '+control-p9_kd-all'
    }

    call macs.macs as SE_nomodel {
        input :
            bamfile=SE_s_mergebam_afterbklist,
            control=SE_c_mergebam_afterbklist,
            nomodel=true,
            egs=egs.genomesize,
            output_name=basename(SE_s_mergebam_afterbklist,'.bam') + '+control-nm',
            default_location = if defined(results_name) then results_name + '/single-end_mode/PEAKS/NARROW_peaks' + '/' + basename(SE_s_mergebam_afterbklist,'.bam') + '+control-nm' else if multi_fastqpair then 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs+control/single-end_mode/PEAKS/NARROW_peaks' + '/' + basename(SE_s_mergebam_afterbklist,'.bam') + '+control-nm' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '+control/single-end_mode/PEAKS/NARROW_peaks' + '/' + basename(SE_s_mergebam_afterbklist,'.bam') + '+control-nm',
            coverage_location = if defined(results_name) then results_name + '/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + basename(SE_s_mergebam_afterbklist,'.bam') + '+control-nm' else if multi_fastqpair then 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs+control/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + basename(SE_s_mergebam_afterbklist,'.bam') + '+control-nm' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '+control/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + basename(SE_s_mergebam_afterbklist,'.bam') + '+control-nm'
    }

    call bamtogff.bamtogff as SE_bamtogff {
        input :
            gtffile=gtf,
            chromsizes=samtools_faidx.chromsizes,
            bamfile=SE_s_merge_markdup.mkdupbam,
            bamindex=SE_s_merge_mkdup.indexbam,
            control_bamfile=SE_c_merge_markdup.mkdupbam,
            control_bamindex=SE_c_merge_mkdup.indexbam,
            samplename=basename(SE_s_merge_markdup.mkdupbam,'.bam') + '+control',
            default_location = if defined(results_name) then results_name + '/single-end_mode/BAM_Density' else if multi_fastqpair then 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs+control/single-end_mode/BAM_Density' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '+control/single-end_mode/BAM_Density'
    }

    call bedtools.bamtobed as SE_s_forsicerbed {
        input :
            bamfile=SE_s_merge_markdup.mkdupbam
    }

    call bedtools.bamtobed as SE_c_forsicerbed {
        input :
            bamfile=SE_c_merge_markdup.mkdupbam
    }

    call sicer.sicer as SE_sicer {
        input :
            bedfile=SE_s_forsicerbed.bedfile,
            control_bed=SE_c_forsicerbed.bedfile,
            chromsizes=samtools_faidx.chromsizes,
            genome_fraction=egs.genomefraction,
            fragmentlength=SE_s_mergebam.avg_readlength,
            outputname=basename(SE_s_forsicerbed.bedfile,'.bed') + '+control',
            default_location = if defined(results_name) then results_name + '/single-end_mode/PEAKS/BROAD_peaks' else if multi_fastqpair then 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs+control/single-end_mode/PEAKS/BROAD_peaks' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '+control/single-end_mode/PEAKS/BROAD_peaks',
            coverage_location = if defined(results_name) then results_name + '/single-end_mode/COVERAGE_files/BROAD_peaks' else if multi_fastqpair then 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs+control/single-end_mode/COVERAGE_files/BROAD_peaks' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '+control/single-end_mode/COVERAGE_files/BROAD_peaks'
    }

    call rose.rose as SE_rose {
        input :
            gtffile=gtf,
            bamfile=SE_s_merge_markdup.mkdupbam,
            bamindex=SE_s_merge_mkdup.indexbam,
            control=SE_c_merge_markdup.mkdupbam,
            controlindex=SE_c_merge_mkdup.indexbam,
            bedfile_auto=SE_macs.peakbedfile,
            bedfile_all=SE_all.peakbedfile,
            default_location = if defined(results_name) then results_name + '/single-end_mode/PEAKS/STITCHED_peaks' else if multi_fastqpair then 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs+control/single-end_mode/PEAKS/STITCHED_peaks' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '+control/single-end_mode/PEAKS/STITCHED_peaks'
    }

    call runspp.runspp as SE_runspp {
        input:
            bamfile=SE_s_mergebam_afterbklist,
            control=SE_c_mergebam_afterbklist
    }

    call runspp.runspp as SE_s_runspp {
        input:
            bamfile=SE_s_mergebam_afterbklist
    }

    call runspp.runspp as SE_c_runspp {
        input:
            bamfile=SE_c_mergebam_afterbklist
    }

    String SE_string_ctrlwig = ""
    call viz.visualization as SE_c_visualization {
        input:
            wigfile=select_first([SE_macs.ctrlwigfile, SE_string_ctrlwig]),
            chromsizes=samtools_faidx.chromsizes,
            control=true,
            xlsfile=SE_macs.peakxlsfile,
            default_location=if defined(results_name) then results_name + '/single-end_mode/COVERAGE_files/NARROW_peaks/' + sub(basename(SE_macs.peakbedfile),'_peaks.bed','') + '/control' else sub(basename(SE_s_mergebam_afterbklist),'.sorted.b.*$','') + '+control/single-end_mode/COVERAGE_files/NARROW_peaks/' + sub(basename(SE_macs.peakbedfile),'_peaks.bed','') + '/control'
    }

    call viz.visualization as SE_c_vizall {
        input:
            wigfile=select_first([SE_all.ctrlwigfile, SE_string_ctrlwig]),
            chromsizes=samtools_faidx.chromsizes,
            control=true,
            xlsfile=SE_all.peakxlsfile,
            default_location=if defined(results_name) then results_name + '/single-end_mode/COVERAGE_files/NARROW_peaks/' + sub(basename(SE_all.peakbedfile),'_peaks.bed','') + '/control' else sub(basename(SE_s_mergebam_afterbklist),'.sorted.b.*$','') + '+control/single-end_mode/COVERAGE_files/NARROW_peaks/' + sub(basename(SE_all.peakbedfile),'_peaks.bed','') + '/control'
    }
    call viz.visualization as SE_c_viznomodel {
        input:
            wigfile=select_first([SE_nomodel.ctrlwigfile, SE_string_ctrlwig]),
            chromsizes=samtools_faidx.chromsizes,
            control=true,
            xlsfile=SE_nomodel.peakxlsfile,
            default_location=if defined(results_name) then results_name + '/single-end_mode/COVERAGE_files/NARROW_peaks/' + sub(basename(SE_nomodel.peakbedfile),'_peaks.bed','') + '/control' else sub(basename(SE_s_mergebam_afterbklist),'.sorted.b.*$','') + '+control/single-end_mode/COVERAGE_files/NARROW_peaks/' + sub(basename(SE_nomodel.peakbedfile),'_peaks.bed','') + '/control'
    }

    call util.peaksanno as SE_peaksanno {
        input :
            gtffile=gtf,
            bedfile=SE_macs.peakbedfile,
            chromsizes=samtools_faidx.chromsizes,
            summitfile=SE_macs.summitsfile,
            default_location = if defined(results_name) then results_name + '/single-end_mode/PEAKS_Annotation/NARROW_peaks' + '/' + sub(basename(SE_macs.peakbedfile),'_peaks.bed','') else if multi_fastqpair then 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs+control/single-end_mode/PEAKS_Annotation/NARROW_peaks' + '/' + sub(basename(SE_macs.peakbedfile),'_peaks.bed','') else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '+control/single-end_mode/PEAKS_Annotation/NARROW_peaks' + '/' + sub(basename(SE_macs.peakbedfile),'_peaks.bed','')
    }

    call util.peaksanno as SE_all_peaksanno {
        input :
            gtffile=gtf,
            bedfile=SE_all.peakbedfile,
            chromsizes=samtools_faidx.chromsizes,
            summitfile=SE_all.summitsfile,
            default_location = if defined(results_name) then results_name + '/single-end_mode/PEAKS_Annotation/NARROW_peaks' + '/' + sub(basename(SE_all.peakbedfile),'_peaks.bed','') else if multi_fastqpair then 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs+control/single-end_mode/PEAKS_Annotation/NARROW_peaks' + '/' + sub(basename(SE_all.peakbedfile),'_peaks.bed','') else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '+control/single-end_mode/PEAKS_Annotation/NARROW_peaks' + '/' + sub(basename(SE_all.peakbedfile),'_peaks.bed','')
    }

    call util.peaksanno as SE_nomodel_peaksanno {
        input :
            gtffile=gtf,
            bedfile=SE_nomodel.peakbedfile,
            chromsizes=samtools_faidx.chromsizes,
            summitfile=SE_nomodel.summitsfile,
            default_location = if defined(results_name) then results_name + '/single-end_mode/PEAKS_Annotation/NARROW_peaks' + '/' + sub(basename(SE_nomodel.peakbedfile),'_peaks.bed','') else if multi_fastqpair then 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs+control/single-end_mode/PEAKS_Annotation/NARROW_peaks' + '/' + sub(basename(SE_nomodel.peakbedfile),'_peaks.bed','') else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '+control/single-end_mode/PEAKS_Annotation/NARROW_peaks' + '/' + sub(basename(SE_nomodel.peakbedfile),'_peaks.bed','')
    }

    call util.peaksanno as SE_sicer_peaksanno {
        input :
            gtffile=gtf,
            bedfile=select_first([SE_sicer.fdrisland, SE_string_ctrlwig]),
            chromsizes=samtools_faidx.chromsizes,
            default_location = if defined(results_name) then results_name + '/single-end_mode/PEAKS_Annotation/BROAD_peaks' else if multi_fastqpair then 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs+control/single-end_mode/PEAKS_Annotation/BROAD_peaks' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '+control/single-end_mode/PEAKS_Annotation/BROAD_peaks'
    }

    # Motif Analysis
    if (run_motifs) {
        call motifs.motifs as SE_motifs {
            input:
                reference=reference,
                reference_index=samtools_faidx.faidx_file,
                bedfile=SE_macs.peakbedfile,
                motif_databases=motif_databases,
                default_location = if defined(results_name) then results_name + '/single-end_mode/MOTIFS' else if multi_fastqpair then 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs+control/single-end_mode/MOTIFS' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '+control/single-end_mode/MOTIFS'
        }

        call util.flankbed as SE_flankbed {
            input :
                bedfile=SE_macs.summitsfile,
                default_location = if defined(results_name) then results_name + '/single-end_mode/MOTIFS' else if multi_fastqpair then 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs+control/single-end_mode/MOTIFS' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '+control/single-end_mode/MOTIFS'
        }

        call motifs.motifs as SE_flank {
            input:
                reference=reference,
                reference_index=samtools_faidx.faidx_file,
                bedfile=SE_flankbed.flankbedfile,
                motif_databases=motif_databases,
                default_location = if defined(results_name) then results_name + '/single-end_mode/MOTIFS' else if multi_fastqpair then 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs+control/single-end_mode/MOTIFS' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '+control/single-end_mode/MOTIFS'
        }
    }

    call viz.visualization as SE_visualization {
        input:
            wigfile=SE_macs.wigfile,
            chromsizes=samtools_faidx.chromsizes,
            xlsfile=SE_macs.peakxlsfile,
            default_location = if defined(results_name) then results_name + '/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + sub(basename(SE_macs.peakbedfile),'_peaks.bed','') else if multi_fastqpair then 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs+control/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + sub(basename(SE_macs.peakbedfile),'_peaks.bed','') else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '+control/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + sub(basename(SE_macs.peakbedfile),'_peaks.bed','')
    }

    call viz.visualization as SE_vizall {
        input:
            wigfile=SE_all.wigfile,
            chromsizes=samtools_faidx.chromsizes,
            xlsfile=SE_all.peakxlsfile,
            default_location = if defined(results_name) then results_name + '/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + sub(basename(SE_all.peakbedfile),'_peaks.bed','') else if multi_fastqpair then 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs+control/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + sub(basename(SE_all.peakbedfile),'_peaks.bed','') else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '+control/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + sub(basename(SE_all.peakbedfile),'_peaks.bed','')
    }

    call viz.visualization as SE_viznomodel {
        input:
            wigfile=SE_nomodel.wigfile,
            chromsizes=samtools_faidx.chromsizes,
            xlsfile=SE_nomodel.peakxlsfile,
            default_location = if defined(results_name) then results_name + '/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + sub(basename(SE_nomodel.peakbedfile),'_peaks.bed','') else if multi_fastqpair then 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs+control/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + sub(basename(SE_nomodel.peakbedfile),'_peaks.bed','') else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '+control/single-end_mode/COVERAGE_files/NARROW_peaks' + '/' + sub(basename(SE_nomodel.peakbedfile),'_peaks.bed','')
    }

    call viz.visualization as SE_vizsicer {
        input:
            wigfile=SE_sicer.wigfile,
            chromsizes=samtools_faidx.chromsizes,
            default_location = if defined(results_name) then results_name + '/single-end_mode/COVERAGE_files/BROAD_peaks' else if multi_fastqpair then 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs+control/single-end_mode/COVERAGE_files/BROAD_peaks' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '+control/single-end_mode/COVERAGE_files/BROAD_peaks'
    }

    #Peak Calling for Sample BAM only
    call macs.macs as only_s_macs {
        input :
            bamfile=SE_s_mergebam_afterbklist,
            pvalue="1e-9",
            keep_dup="auto",
            egs=egs.genomesize,
            default_location=if multi_fastqpair then 'SAMPLE/AllCases_' + length(sample_fastqfiles) + 'fastqpairs/single-end_mode/PEAKS_forQC/' + basename(SE_s_mergebam_afterbklist,'.bam') + '-p9_kd-auto' else 'SAMPLE/' + sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/PEAKS_forQC/' + basename(SE_s_mergebam_afterbklist,'.bam') + '-p9_kd-auto',
            coverage_location=if multi_fastqpair then 'SAMPLE/AllCases_' + length(sample_fastqfiles) + 'fastqpairs/single-end_mode/PEAKS_forQC/' + basename(SE_s_mergebam_afterbklist,'.bam') + '-p9_kd-auto' else 'SAMPLE/' + sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/PEAKS_forQC/' + basename(SE_s_mergebam_afterbklist,'.bam') + '-p9_kd-auto'

    }

    #Peak Calling for Control BAM only
    call macs.macs as only_c_macs {
        input :
            bamfile=SE_c_mergebam_afterbklist,
            pvalue="1e-9",
            keep_dup="auto",
            egs=egs.genomesize,
            default_location=if multi_control_fastqpair then 'CONTROL/AllControls_' + length(control_fastqfiles) + 'fastqpairs/single-end_mode/PEAKS_forQC/' + basename(SE_c_mergebam_afterbklist,'.bam') + '-p9_kd-auto' else 'CONTROL/' + sub(basename(control_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/PEAKS_forQC/' + basename(SE_c_mergebam_afterbklist,'.bam') + '-p9_kd-auto',
            coverage_location=if multi_control_fastqpair then 'CONTROL/AllControls_' + length(control_fastqfiles) + 'fastqpairs/single-end_mode/PEAKS_forQC/' + basename(SE_c_mergebam_afterbklist,'.bam') + '-p9_kd-auto' else 'CONTROL/' + sub(basename(control_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/PEAKS_forQC/' + basename(SE_c_mergebam_afterbklist,'.bam') + '-p9_kd-auto'

    }


    call bedtools.bamtobed as SE_c_finalbed {
        input:
            bamfile=SE_c_mergebam_afterbklist
    }

    call sortbed.sortbed as SE_c_sortbed {
        input:
            bedfile=SE_c_finalbed.bedfile
    }

    call bedtools.intersect as SE_c_intersect {
        input:
            fileA=only_c_macs.peakbedfile,
            fileB=SE_c_sortbed.sortbed_out,
            countoverlap=true,
            sorted=true
    }

    call bedtools.bamtobed as SE_s_finalbed {
        input:
            bamfile=SE_s_mergebam_afterbklist
    }

    call sortbed.sortbed as SE_s_sortbed {
        input:
            bedfile=SE_s_finalbed.bedfile
    }

    call bedtools.intersect as SE_s_intersect {
        input:
            fileA=only_s_macs.peakbedfile,
            fileB=SE_s_sortbed.sortbed_out,
            countoverlap=true,
            sorted=true
    }

    call bedtools.intersect as SE_intersect {
        input:
            fileA=SE_macs.peakbedfile,
            fileB=SE_s_sortbed.sortbed_out,
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
    File PE_sample_bam = select_first([s_PE_mergebam_afterbklist, s_uno_PE_mapping.bklist_bam, s_uno_PE_mapping.sorted_bam])
    File PE_control_bam = select_first([c_PE_mergebam_afterbklist, c_uno_PE_mapping.bklist_bam, c_uno_PE_mapping.sorted_bam])

    call macs.macs as PE_macs {
        input :
            bamfile=PE_sample_bam,
            control=PE_control_bam,
            pvalue="1e-9",
            keep_dup="auto",
            egs=egs.genomesize,
            output_name = if defined(results_name) then results_name + '-p9_kd-auto' else basename(PE_sample_bam,'.bam') + '+control-p9_kd-auto',
            default_location=if defined(results_name) then results_name + '/PEAKS/NARROW_peaks/' + results_name + '-p9_kd-auto' else sub(basename(PE_sample_bam),'.sorted.b.*$','') + '+control/PEAKS/NARROW_peaks' + '/' + basename(PE_sample_bam,'.bam') + '+control-p9_kd-auto',
            coverage_location=if defined(results_name) then results_name + '/COVERAGE_files/NARROW_peaks/' + results_name + '-p9_kd-auto' else sub(basename(PE_sample_bam),'.sorted.b.*$','') + '+control/COVERAGE_files/NARROW_peaks' + '/' + basename(PE_sample_bam,'.bam') + '+control-p9_kd-auto'
    }

    call util.addreadme as PE_addreadme {
        input :
            default_location=if defined(results_name) then results_name + '/PEAKS' else sub(basename(PE_sample_bam),'.sorted.b.*$','') + '+control/PEAKS'
    }

    call macs.macs as PE_all {
        input :
            bamfile=PE_sample_bam,
            control=PE_control_bam,
            pvalue="1e-9",
            keep_dup="all",
            egs=egs.genomesize,
            output_name = if defined(results_name) then results_name + '-p9_kd-all' else basename(PE_sample_bam,'.bam') + '+control-p9_kd-all',
            default_location=if defined(results_name) then results_name + '/PEAKS/NARROW_peaks/' + results_name + '-p9_kd-all' else sub(basename(PE_sample_bam),'.sorted.b.*$','') + '+control/PEAKS/NARROW_peaks' + '/' + basename(PE_sample_bam,'.bam') + '+control-p9_kd-all',
            coverage_location=if defined(results_name) then results_name + '/COVERAGE_files/NARROW_peaks/' + results_name + '-p9_kd-all' else sub(basename(PE_sample_bam),'.sorted.b.*$','') + '+control/COVERAGE_files/NARROW_peaks' + '/' + basename(PE_sample_bam,'.bam') + '+control-p9_kd-all'
    }

    call macs.macs as PE_nomodel {
        input :
            bamfile=PE_sample_bam,
            control=PE_control_bam,
            nomodel=true,
            egs=egs.genomesize,
            output_name=if defined(results_name) then results_name + '-nm' else basename(PE_sample_bam,'.bam') + '+control-nm',
            default_location=if defined(results_name) then results_name + '/PEAKS/NARROW_peaks/' + results_name + '-nm' else sub(basename(PE_sample_bam),'.sorted.b.*$','') + '+control/PEAKS/NARROW_peaks' + '/' + basename(PE_sample_bam,'.bam') + '+control-nm',
            coverage_location=if defined(results_name) then results_name + '/COVERAGE_files/NARROW_peaks/' + results_name + '-nm' else sub(basename(PE_sample_bam),'.sorted.b.*$','') + '+control/COVERAGE_files/NARROW_peaks' + '/' + basename(PE_sample_bam,'.bam') + '+control-nm'
    }

    call peaseq_util.fraggraph as s_fraggraph {
        input :
            bamfile=select_first([s_PE_merge_markdup.mkdupbam, s_uno_PE_mapping.mkdup_bam]),
            chromsizes=samtools_faidx.chromsizes,
            default_location='SAMPLE/' + sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/COVERAGE_files/BAM_Fragments',
            bam_location='SAMPLE/' + sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/BAM_files',
            annotation_location='SAMPLE/' + sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/QC/Fragments'
    }

    call samtools.indexstats as s_frag_index {
        input :
            bamfile=s_fraggraph.fragbamfile,
            default_location='SAMPLE/' + sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/BAM_files'
    }

    call peaseq_util.fraggraph as c_fraggraph {
        input :
            bamfile=select_first([c_PE_merge_markdup.mkdupbam, c_uno_PE_mapping.mkdup_bam]),
            chromsizes=samtools_faidx.chromsizes,
            default_location='CONTROL/' + sub(basename(PE_control_bam),'.sorted.b.*$','') + '/COVERAGE_files/BAM_Fragments',
            bam_location='CONTROL/' + sub(basename(PE_control_bam),'.sorted.b.*$','') + '/BAM_files',
            annotation_location='CONTROL/' + sub(basename(PE_control_bam),'.sorted.b.*$','') + '/QC/Fragments'
    }

    call samtools.indexstats as c_frag_index {
        input :
            bamfile=c_fraggraph.fragbamfile,
            default_location='CONTROL/' + sub(basename(PE_control_bam),'.sorted.b.*$','') + '/BAM_files'
    }

    call bamtogff.bamtogff as PE_bamtogff {
        input :
            gtffile=gtf,
            chromsizes=samtools_faidx.chromsizes,
            bamfile=s_fraggraph.fragbamfile,
            bamindex=s_frag_index.indexbam,
            control_bamfile=c_fraggraph.fragbamfile,
            control_bamindex=c_frag_index.indexbam,
            samplename=if defined(results_name) then results_name else basename(PE_sample_bam,'.bam') + '+control',
            default_location=if defined(results_name) then results_name + '/BAM_Density' else sub(basename(PE_sample_bam),'.sorted.b.*$','') + '+control/BAM_Density'
    }

    call sicer.sicer as PE_sicer {
        input :
            bedfile=s_fraggraph.bedpefile,
            control_bed=c_fraggraph.bedpefile,
            paired_end=true,
            gap_size=600,
            chromsizes=samtools_faidx.chromsizes,
            genome_fraction=egs.genomefraction,
            outputname=if defined(results_name) then results_name else basename(s_fraggraph.bedpefile,'.bed') + '+control',
            default_location=if defined(results_name) then results_name + '/PEAKS/BROAD_peaks' else sub(basename(PE_sample_bam),'.sorted.b.*$','') + '+control/PEAKS/BROAD_peaks',
            coverage_location=if defined(results_name) then results_name + '/PEAKS/BROAD_peaks' else sub(basename(PE_sample_bam),'.sorted.b.*$','') + '+control/COVERAGE_files/BROAD_peaks'
    }

    call rose.rose as PE_rose {
        input :
            gtffile=gtf,
            bamfile=s_fraggraph.fragbamfile,
            bamindex=s_frag_index.indexbam,
            control=c_fraggraph.fragbamfile,
            controlindex=c_frag_index.indexbam,
            bedfile_auto=PE_macs.peakbedfile,
            bedfile_all=PE_all.peakbedfile,
            default_location=if defined(results_name) then results_name + '/PEAKS/STITCHED_peaks' else sub(basename(PE_sample_bam),'.sorted.b.*$','') + '+control/PEAKS/STITCHED_peaks'
    }

    call runspp.runspp as PE_runspp {
        input:
            bamfile=PE_sample_bam,
            control=PE_control_bam
    }

    call runspp.runspp as s_PE_runspp {
        input:
            bamfile=PE_sample_bam
    }

    call runspp.runspp as c_PE_runspp {
        input:
            bamfile=PE_control_bam
    }

    String string_ctrlwig = ""
    call viz.visualization as c_PE_visualization {
        input:
            wigfile=select_first([PE_macs.ctrlwigfile, string_ctrlwig]),
            chromsizes=samtools_faidx.chromsizes,
            control=true,
            xlsfile=PE_macs.peakxlsfile,
            default_location=if defined(results_name) then results_name + '/COVERAGE_files/NARROW_peaks/' + sub(basename(PE_macs.peakbedfile),'_peaks.bed','') + '/control' else sub(basename(PE_sample_bam),'.sorted.b.*$','') + '+control/COVERAGE_files/NARROW_peaks/' + sub(basename(PE_macs.peakbedfile),'_peaks.bed','') + '/control'
    }

    call viz.visualization as c_PE_vizall {
        input:
            wigfile=select_first([PE_all.ctrlwigfile, string_ctrlwig]),
            chromsizes=samtools_faidx.chromsizes,
            control=true,
            xlsfile=PE_all.peakxlsfile,
            default_location=if defined(results_name) then results_name + '/COVERAGE_files/NARROW_peaks/' + sub(basename(PE_all.peakbedfile),'_peaks.bed','') + '/control' else sub(basename(PE_sample_bam),'.sorted.b.*$','') + '+control/COVERAGE_files/NARROW_peaks/' + sub(basename(PE_all.peakbedfile),'_peaks.bed','') + '/control'
    }
    call viz.visualization as c_PE_viznomodel {
        input:
            wigfile=select_first([PE_nomodel.ctrlwigfile, string_ctrlwig]),
            chromsizes=samtools_faidx.chromsizes,
            control=true,
            xlsfile=PE_nomodel.peakxlsfile,
            default_location=if defined(results_name) then results_name + '/COVERAGE_files/NARROW_peaks/' + sub(basename(PE_nomodel.peakbedfile),'_peaks.bed','') + '/control' else sub(basename(PE_sample_bam),'.sorted.b.*$','') + '+control/COVERAGE_files/NARROW_peaks/' + sub(basename(PE_nomodel.peakbedfile),'_peaks.bed','') + '/control'
    }

    call util.peaksanno as PE_peaksanno {
        input :
            gtffile=gtf,
            bedfile=PE_macs.peakbedfile,
            chromsizes=samtools_faidx.chromsizes,
            summitfile=PE_macs.summitsfile,
            default_location=if defined(results_name) then results_name + '/PEAKS_Annotation/NARROW_peaks/' + sub(basename(PE_macs.peakbedfile),'_peaks.bed','') else sub(basename(PE_sample_bam),'.sorted.b.*$','') + '+control/PEAKS_Annotation/NARROW_peaks/' + sub(basename(PE_macs.peakbedfile),'_peaks.bed','')
    }

    call util.peaksanno as PE_all_peaksanno {
        input :
            gtffile=gtf,
            bedfile=PE_all.peakbedfile,
            chromsizes=samtools_faidx.chromsizes,
            summitfile=PE_all.summitsfile,
            default_location=if defined(results_name) then results_name + '/PEAKS_Annotation/NARROW_peaks/' + sub(basename(PE_all.peakbedfile),'_peaks.bed','') else sub(basename(PE_sample_bam),'.sorted.b.*$','') + '+control/PEAKS_Annotation/NARROW_peaks/' + sub(basename(PE_all.peakbedfile),'_peaks.bed','')
    }

    call util.peaksanno as PE_nomodel_peaksanno {
        input :
            gtffile=gtf,
            bedfile=PE_nomodel.peakbedfile,
            chromsizes=samtools_faidx.chromsizes,
            summitfile=PE_nomodel.summitsfile,
            default_location=if defined(results_name) then results_name + '/PEAKS_Annotation/NARROW_peaks/' + sub(basename(PE_nomodel.peakbedfile),'_peaks.bed','') else sub(basename(PE_sample_bam),'.sorted.b.*$','') + '+control/PEAKS_Annotation/NARROW_peaks/' + sub(basename(PE_nomodel.peakbedfile),'_peaks.bed','')
    }

    call util.peaksanno as PE_sicer_peaksanno {
        input :
            gtffile=gtf,
            bedfile=select_first([PE_sicer.fdrisland, string_ctrlwig]),
            chromsizes=samtools_faidx.chromsizes,
            default_location=if defined(results_name) then results_name + '/PEAKS_Annotation/BROAD_peaks' else sub(basename(PE_sample_bam),'.sorted.b.*$','') + '+control/PEAKS_Annotation/BROAD_peaks'
    }

    # Motif Analysis
    if (run_motifs) {
        call motifs.motifs as PE_motifs {
            input:
                reference=reference,
                reference_index=samtools_faidx.faidx_file,
                bedfile=PE_macs.peakbedfile,
                motif_databases=motif_databases,
                default_location=if defined(results_name) then results_name + '/MOTIFS' else sub(basename(PE_sample_bam),'.sorted.b.*$','') + '+control/MOTIFS'
        }

        call util.flankbed as PE_flankbed {
            input :
                bedfile=PE_macs.summitsfile,
                default_location=if defined(results_name) then results_name + '/MOTIFS' else sub(basename(PE_sample_bam),'.sorted.b.*$','') + '+control/MOTIFS'
        }

        call motifs.motifs as PE_flank {
            input:
                reference=reference,
                reference_index=samtools_faidx.faidx_file,
                bedfile=PE_flankbed.flankbedfile,
                motif_databases=motif_databases,
                default_location=if defined(results_name) then results_name + '/MOTIFS' else sub(basename(PE_sample_bam),'.sorted.b.*$','') + '+control/MOTIFS'
        }
    }

    call viz.visualization as PE_visualization {
        input:
            wigfile=PE_macs.wigfile,
            chromsizes=samtools_faidx.chromsizes,
            xlsfile=PE_macs.peakxlsfile,
            default_location=if defined(results_name) then results_name + '/COVERAGE_files/NARROW_peaks/' + sub(basename(PE_macs.peakbedfile),'_peaks.bed','') else sub(basename(PE_sample_bam),'.sorted.b.*$','') + '+control/COVERAGE_files/NARROW_peaks/' + sub(basename(PE_macs.peakbedfile),'_peaks.bed','')
    }

    call viz.visualization as PE_vizall {
        input:
            wigfile=PE_all.wigfile,
            chromsizes=samtools_faidx.chromsizes,
            xlsfile=PE_all.peakxlsfile,
            default_location=if defined(results_name) then results_name + '/COVERAGE_files/NARROW_peaks/' + sub(basename(PE_all.peakbedfile),'_peaks.bed','') else sub(basename(PE_sample_bam),'.sorted.b.*$','') + '+control/COVERAGE_files/NARROW_peaks/' + sub(basename(PE_all.peakbedfile),'_peaks.bed','')
    }

    call viz.visualization as PE_viznomodel {
        input:
            wigfile=PE_nomodel.wigfile,
            chromsizes=samtools_faidx.chromsizes,
            xlsfile=PE_nomodel.peakxlsfile,
            default_location=if defined(results_name) then results_name + '/COVERAGE_files/NARROW_peaks/' + sub(basename(PE_nomodel.peakbedfile),'_peaks.bed','') else sub(basename(PE_sample_bam),'.sorted.b.*$','') + '+control/COVERAGE_files/NARROW_peaks/' + sub(basename(PE_nomodel.peakbedfile),'_peaks.bed','')
    }

    call viz.visualization as PE_vizsicer {
        input:
            wigfile=PE_sicer.wigfile,
            chromsizes=samtools_faidx.chromsizes,
            default_location=if defined(results_name) then results_name + '/COVERAGE_files/BROAD_peaks' else sub(basename(PE_sample_bam),'.sorted.b.*$','') + '+control/COVERAGE_files/BROAD_peaks'
    }

    #Peak Calling for Sample BAM only
    call macs.macs as only_s_PE_macs {
        input :
            bamfile=PE_sample_bam,
            pvalue="1e-9",
            keep_dup="auto",
            egs=egs.genomesize,
            default_location='SAMPLE/' + sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/PEAKS_forQC/' + basename(PE_sample_bam,'.bam') + '-p9_kd-auto',
            coverage_location='SAMPLE/' + sub(basename(PE_sample_bam),'.sorted.b.*$','') + '/PEAKS_forQC/' + basename(PE_sample_bam,'.bam') + '-p9_kd-auto'
    }

    #Peak Calling for Control BAM only
    call macs.macs as only_c_PE_macs {
        input :
            bamfile=PE_control_bam,
            pvalue="1e-9",
            keep_dup="auto",
            egs=egs.genomesize,
            default_location='CONTROL/' + sub(basename(PE_control_bam),'.sorted.b.*$','') + '/PEAKS_forQC/' + basename(PE_control_bam,'.bam') + '-p9_kd-auto',
            coverage_location='CONTROL/' + sub(basename(PE_control_bam),'.sorted.b.*$','') + '/PEAKS_forQC/' + basename(PE_control_bam,'.bam') + '-p9_kd-auto'
    }

    call bedtools.bamtobed as only_c_PE_finalbed {
        input:
            bamfile=PE_control_bam
    }

    call sortbed.sortbed as only_c_PE_sortbed {
        input:
            bedfile=only_c_PE_finalbed.bedfile
    }

    call bedtools.intersect as only_c_PE_intersect {
        input:
            fileA=only_c_PE_macs.peakbedfile,
            fileB=only_c_PE_sortbed.sortbed_out,
            countoverlap=true,
            sorted=true
    }

    call bedtools.bamtobed as only_s_PE_finalbed {
        input:
            bamfile=PE_sample_bam
    }

    call sortbed.sortbed as only_s_PE_sortbed {
        input:
            bedfile=only_s_PE_finalbed.bedfile
    }

    call bedtools.intersect as only_s_PE_intersect {
        input:
            fileA=only_s_PE_macs.peakbedfile,
            fileB=only_s_PE_sortbed.sortbed_out,
            countoverlap=true,
            sorted=true
    }

    call bedtools.intersect as PE_intersect {
        input:
            fileA=PE_macs.peakbedfile,
            fileB=only_s_PE_sortbed.sortbed_out,
            countoverlap=true,
            sorted=true
    }

    # Final MERGE output file
    File s_mergehtmlfile =  select_first([s_final_mergehtml.mergefile, SE_s_mergehtml.mergefile])
    File c_mergehtmlfile =  select_first([c_final_mergehtml.mergefile, SE_c_mergehtml.mergefile])

### ------------------------------------------------- ###
### ----------------- S E C T I O N 5 --------------- ###
### --------------- Summary Statistics -------------- ###
### ------------ A: analysis for SE mode  ----------- ###
### ------------------------------------------------- ###
    call util.evalstats as SE_s_summarystats {
        input:
            fastq_type="PEAseq SEmode Sample",
            bambed=SE_s_finalbed.bedfile,
            sppfile=SE_s_runspp.spp_out,
            fastqczip=SE_s_mergebamfqc.zipfile,
            bamflag=SE_s_mergeindexstats.flagstats,
            rmdupflag=SE_s_merge_mkdup.flagstats,
            bkflag=SE_s_merge_bklist.flagstats,
            countsfile=SE_s_intersect.intersect_out,
            peaksxls=only_s_macs.peakxlsfile,
            outputfile = sub(basename(SE_s_mergebam_afterbklist),'.sorted.b.*$', '-stats.csv'),
            outputhtml = sub(basename(SE_s_mergebam_afterbklist),'.sorted.b.*$','-stats.html'),
            outputtext = sub(basename(SE_s_mergebam_afterbklist),'.sorted.b.*$', '-stats.txt'),
            configml = sub(basename(SE_s_mergebam_afterbklist),'.sorted.b.*$', '-config.ml'),
            default_location = if multi_fastqpair then 'SAMPLE/AllCases_' + length(sample_fastqfiles) + 'fastqpairs/single-end_mode/QC/SummaryStats' else 'SAMPLE/' + sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/QC/SummaryStats'
    }

    call util.evalstats as merge_SE_summarystats {
        input:
            bambed=SE_s_finalbed.bedfile,
            sppfile=SE_runspp.spp_out,
            fastqczip=SE_s_mergebamfqc.zipfile,
            bamflag=SE_s_mergeindexstats.flagstats,
            rmdupflag=SE_s_merge_mkdup.flagstats,
            bkflag=SE_s_merge_bklist.flagstats,
            countsfile=SE_intersect.intersect_out,
            peaksxls=SE_macs.peakxlsfile,
            enhancers=SE_rose.enhancers,
            superenhancers=SE_rose.super_enhancers,
            outputfile = sub(basename(SE_s_mergebam_afterbklist),'.sorted.b.*$','-stats.csv'),
            outputhtml = sub(basename(SE_s_mergebam_afterbklist),'.sorted.b.*$','-stats.html'),
            outputtext = sub(basename(SE_s_mergebam_afterbklist),'.sorted.b.*$','-stats.txt'),
            configml = sub(basename(SE_s_mergebam_afterbklist),'.sorted.b.*$','-config.ml')
    }


    call util.evalstats as SE_c_summarystats {
        input:
            fastq_type="PEAseq SEmode Control",
            bambed=SE_c_finalbed.bedfile,
            sppfile=SE_c_runspp.spp_out,
            fastqczip=SE_c_mergebamfqc.zipfile,
            bamflag=SE_c_mergeindexstats.flagstats,
            rmdupflag=SE_c_merge_mkdup.flagstats,
            bkflag=SE_c_merge_bklist.flagstats,
            countsfile=SE_c_intersect.intersect_out,
            peaksxls=only_c_macs.peakxlsfile,
            outputfile = sub(basename(SE_c_mergebam_afterbklist),'.sorted.b.*$', '-stats.csv'),
            outputhtml = sub(basename(SE_c_mergebam_afterbklist),'.sorted.b.*$','-stats.html'),
            outputtext = sub(basename(SE_c_mergebam_afterbklist),'.sorted.b.*$', '-stats.txt'),
            configml = sub(basename(SE_c_mergebam_afterbklist),'.sorted.b.*$', '-config.ml'),
            default_location = if multi_control_fastqpair then 'CONTROL/AllControls_' + length(control_fastqfiles) + 'fastqpairs/single-end_mode/QC/SummaryStats' else 'CONTROL/' + sub(basename(control_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/single-end_mode/QC/SummaryStats'
    }

    call util.concatstats as SE_concatstats {
        # SUMMARY STATISTICS of all samples files (more than 1 sample file provided)
        input:
            fastq_mode="PEAseq SEmode",
            sample_config=SE_s_summarystats.configfile,
            control_config=SE_c_summarystats.configfile,
            overall_config=merge_SE_summarystats.configfile,
            outputfile='AllCases_' + length(all_sample_fastqfiles) + 'fastqs+control',
            default_location=if defined(results_name) then results_name + '/single-end_mode/QC/SummaryStats' else if multi_fastqpair then 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs+control/single-end_mode/QC/SummaryStats' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '+control/single-end_mode/QC/SummaryStats'
    }

    call util.summaryreport as SE_overallsummary {
        # Presenting all quality stats for the analysis
        input:
            controlqc_html=SE_c_mergehtml.xhtml,
            sampleqc_html=SE_s_mergehtml.xhtml,
            overallqc_html=SE_concatstats.xhtml,
            controlqc_txt=SE_c_mergehtml.mergetxt,
            sampleqc_txt=SE_s_mergehtml.mergetxt,
            overallqc_txt=SE_concatstats.textfile,
            fastq_mode="PEAseq SEmode",
            default_location=if defined(results_name) then results_name + '/single-end_mode' else if multi_fastqpair then 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs+control/single-end_mode' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '+control/single-end_mode'
    }

### ------------------------------------------------- ###
### ----------------- S E C T I O N 5 --------------- ###
### --------------- Summary Statistics -------------- ###
### ------------ B: analysis for PE mode  ----------- ###
### ------------------------------------------------- ###

    String string_qual = "" #buffer to allow for optionality in if statement

    if ( one_fastqpair ) {
        call util.evalstats as PE_s_summarystats {
            input:
                fastq_type="PEAseq PEmode Sample",
                bambed=only_s_PE_finalbed.bedfile,
                sppfile=s_PE_runspp.spp_out,
                fastqczip=select_first([s_uno_PE_bamfqc.zipfile, string_qual]),
                bamflag=s_uno_PE_mapping.bam_stats,
                rmdupflag=s_uno_PE_mapping.mkdup_stats,
                bkflag=s_uno_PE_mapping.bklist_stats,
                fastqmetrics=indv_s_bfs.metrics_out[0],
                countsfile=only_s_PE_intersect.intersect_out,
                peaksxls=only_s_PE_macs.peakxlsfile,
                outputfile = sub(basename(PE_sample_bam),'.sorted.b.*$', '-stats.csv'),
                outputhtml = sub(basename(PE_sample_bam),'.sorted.b.*$','-stats.html'),
                outputtext = sub(basename(PE_sample_bam),'.sorted.b.*$', '-stats.txt'),
                configml = sub(basename(PE_sample_bam),'.sorted.b.*$', '-config.ml'),
                default_location = 'SAMPLE/' + sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/SummaryStats'
        }

        call util.evalstats as PE_all_summarystats {
            # SUMMARY STATISTICS of sample file (only 1 sample file provided)
            input:
                fastq_type="PEAseq PEmode Comprehensive",
                bambed=only_s_PE_finalbed.bedfile,
                sppfile=PE_runspp.spp_out,
                fastqczip=select_first([s_uno_PE_bamfqc.zipfile, string_qual]),
                bamflag=s_uno_PE_mapping.bam_stats,
                rmdupflag=s_uno_PE_mapping.mkdup_stats,
                bkflag=s_uno_PE_mapping.bklist_stats,
                fastqmetrics=indv_s_bfs.metrics_out[0],
                countsfile=PE_intersect.intersect_out,
                peaksxls=PE_macs.peakxlsfile,
                enhancers=PE_rose.enhancers,
                superenhancers=PE_rose.super_enhancers,
                outputfile = sub(basename(PE_sample_bam),'.sorted.b.*$', '-stats.csv'),
                outputhtml = sub(basename(PE_sample_bam),'.sorted.b.*$','-stats.html'),
                outputtext = sub(basename(PE_sample_bam),'.sorted.b.*$', '-stats.txt'),
                configml = sub(basename(PE_sample_bam),'.sorted.b.*$', '-config.ml')
        }
    } #if one_fastqpair

    if (one_control_fastqpair) {
        call util.evalstats as PE_c_summarystats {
            input:
                fastq_type="PEAseq PEmode Control",
                bambed=only_c_PE_finalbed.bedfile,
                sppfile=c_PE_runspp.spp_out,
                fastqczip=select_first([c_uno_PE_bamfqc.zipfile, string_qual]),
                bamflag=c_uno_PE_mapping.bam_stats,
                rmdupflag=c_uno_PE_mapping.mkdup_stats,
                bkflag=c_uno_PE_mapping.bklist_stats,
                fastqmetrics=indv_c_bfs.metrics_out[0],
                countsfile=only_c_PE_intersect.intersect_out,
                peaksxls=only_c_PE_macs.peakxlsfile,
                outputfile = sub(basename(PE_control_bam),'.sorted.b.*$', '-stats.csv'),
                outputhtml = sub(basename(PE_control_bam),'.sorted.b.*$','-stats.html'),
                outputtext = sub(basename(PE_control_bam),'.sorted.b.*$', '-stats.txt'),
                configml = sub(basename(PE_control_bam),'.sorted.b.*$', '-config.ml'),
                default_location = 'CONTROL/' + sub(basename(control_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '/QC/SummaryStats'
        }
    } #if one_fastqpair

    if ( multi_fastqpair ) {
        call util.evalstats as PE_merge_s_summarystats {
            # SUMMARY STATISTICS of all samples files (more than 1 sample file provided)
            input:
                fastq_type="PEAseq PEmode Sample",
                bambed=only_s_PE_finalbed.bedfile,
                sppfile=s_PE_runspp.spp_out,
                fastqczip=select_first([s_PE_mergebamfqc.zipfile, string_qual]),
                bamflag=s_PE_mergeindexstats.flagstats,
                rmdupflag=s_PE_merge_mkdup.flagstats,
                bkflag=s_PE_merge_bklist.flagstats,
                countsfile=only_s_PE_intersect.intersect_out,
                peaksxls=only_s_PE_macs.peakxlsfile,
                outputfile = sub(basename(PE_sample_bam),'.sorted.b.*$', '-stats.csv'),
                outputhtml = sub(basename(PE_sample_bam),'.sorted.b.*$','-stats.html'),
                outputtext = sub(basename(PE_sample_bam),'.sorted.b.*$', '-stats.txt'),
                configml = sub(basename(PE_sample_bam),'.sorted.b.*$', '-config.ml'),
                default_location = 'SAMPLE/AllCases_' + length(sample_fastqfiles) + 'fastqpairs/QC/SummaryStats'
        }

        call util.evalstats as PE_merge_summarystats {
            # SUMMARY STATISTICS of sample file (only 1 sample file provided)
            input:
                fastq_type="PEAseq PEmode Comprehensive",
                bambed=only_s_PE_finalbed.bedfile,
                sppfile=PE_runspp.spp_out,
                fastqczip=select_first([s_PE_mergebamfqc.zipfile, string_qual]),
                bamflag=s_PE_mergeindexstats.flagstats,
                rmdupflag=s_PE_merge_mkdup.flagstats,
                bkflag=s_PE_merge_bklist.flagstats,
                countsfile=PE_intersect.intersect_out,
                peaksxls=PE_macs.peakxlsfile,
                enhancers=PE_rose.enhancers,
                superenhancers=PE_rose.super_enhancers,
                outputfile = sub(basename(PE_sample_bam),'.sorted.b.*$', '-stats.csv'),
                outputhtml = sub(basename(PE_sample_bam),'.sorted.b.*$','-stats.html'),
                outputtext = sub(basename(PE_sample_bam),'.sorted.b.*$', '-stats.txt'),
                configml = sub(basename(PE_sample_bam),'.sorted.b.*$', '-config.ml')
        }
    } # end if multifastqpair

    if ( multi_control_fastqpair ) {
        call util.evalstats as PE_merge_c_summarystats {
            input:
                fastq_type="PEAseq PEmode Control",
                bambed=only_c_PE_finalbed.bedfile,
                sppfile=c_PE_runspp.spp_out,
                fastqczip=select_first([c_PE_mergebamfqc.zipfile, string_qual]),
                bamflag=c_PE_mergeindexstats.flagstats,
                rmdupflag=c_PE_merge_mkdup.flagstats,
                bkflag=c_PE_merge_bklist.flagstats,
                countsfile=only_c_PE_intersect.intersect_out,
                peaksxls=only_c_PE_macs.peakxlsfile,
                outputfile = sub(basename(PE_control_bam),'.sorted.b.*$', '-stats.csv'),
                outputhtml = sub(basename(PE_control_bam),'.sorted.b.*$','-stats.html'),
                outputtext = sub(basename(PE_control_bam),'.sorted.b.*$', '-stats.txt'),
                configml = sub(basename(PE_control_bam),'.sorted.b.*$', '-config.ml'),
                default_location = 'CONTROL/AllControls_' + length(control_fastqfiles) + 'fastqpairs/QC/SummaryStats'
        }
    } #if one_fastqpair

    call util.concatstats as PE_concatstats {
        # SUMMARY STATISTICS of all samples files (more than 1 sample file provided)
        input:
            fastq_mode="PEAseq PEmode",
            peaseq=true,
            sample_config=select_first([PE_s_summarystats.configfile, PE_merge_s_summarystats.configfile]),
            control_config=select_first([PE_c_summarystats.configfile, PE_merge_c_summarystats.configfile]),
            overall_config=select_first([PE_all_summarystats.configfile, PE_merge_summarystats.configfile]),
            outputfile=if defined(results_name) then results_name else 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs+control',
            default_location=if defined(results_name) then results_name + '/QC/SummaryStats' else if multi_fastqpair then 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs+control/QC/SummaryStats' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '+control/QC/SummaryStats'
    }

    call peaseq_util.pairedend_summaryreport as PE_overallsummary {
        # Presenting all quality stats for the analysis
        input:
            sampleqc_se_html=SE_s_mergehtml.xhtml,
            controlqc_se_html=SE_c_mergehtml.xhtml,
            overallqc_se_html=SE_concatstats.xhtml,
            sampleqc_se_txt=SE_s_mergehtml.mergetxt,
            controlqc_se_txt=SE_c_mergehtml.mergetxt,
            overallqc_se_txt=SE_concatstats.textfile,
            sampleqc_pe_html=select_first([s_PE_mergehtml.xhtml,PE_s_summarystats.xhtml]),
            controlqc_pe_html=select_first([c_PE_mergehtml.xhtml,PE_c_summarystats.xhtml]),
            overallqc_pe_html=PE_concatstats.xhtml,
            sampleqc_pe_txt=select_first([s_PE_mergehtml.mergetxt,PE_s_summarystats.textfile]),
            controlqc_pe_txt=select_first([c_PE_mergehtml.mergetxt,PE_c_summarystats.textfile]),
            overallqc_pe_txt=PE_concatstats.textfile,
            outputfile = if defined(results_name) then results_name + '.peaseq_report.html' else if multi_fastqpair then 'AllCases_' + length(sample_fastqfiles) + 'fastqpairs+control.peaseq_report.html' else sub(basename(sample_fastqfiles[0].left),'_R?[12]_....f.*q.gz|_R?[12].f.*q.gz','') + '+control.peaseq_report.html'
    }

### ------------------------------------------------- ###
### ---------------- S E C T I O N 6 ---------------- ###
### ------------------ OUTPUT FILES ----------------- ###
### ------------------------------------------------- ###

    output {
        #SPIKE-IN
        Array[File?]? spikein_indv_s_htmlfile = spikein_s_indv_fastqc.htmlfile
        Array[File?]? spikein_indv_s_zipfile = spikein_s_indv_fastqc.zipfile
        Array[File?]? spikein_s_metrics_out = spikein_s_indv_map.mapping_output
        Array[File?]? spikein_indv_c_htmlfile = spikein_c_indv_fastqc.htmlfile
        Array[File?]? spikein_indv_c_zipfile = spikein_c_indv_fastqc.zipfile
        Array[File?]? spikein_c_metrics_out = spikein_c_indv_map.mapping_output

        #FASTQC
        Array[File?]? indv_s_htmlfile = s_indv_fastqc.htmlfile
        Array[File?]? indv_s_zipfile = s_indv_fastqc.zipfile
        Array[File?]? indv_s_bam_htmlfile = indv_s_bamfqc.htmlfile
        Array[File?]? indv_s_bam_zipfile = indv_s_bamfqc.zipfile
        File? s_mergebam_htmlfile = SE_s_mergebamfqc.htmlfile
        File? s_mergebam_zipfile = SE_s_mergebamfqc.zipfile

        Array[File?]? indv_c_htmlfile = c_indv_fastqc.htmlfile
        Array[File?]? indv_c_zipfile = c_indv_fastqc.zipfile
        Array[File?]? indv_c_bam_htmlfile = indv_c_bamfqc.htmlfile
        Array[File?]? indv_c_bam_zipfile = indv_c_bamfqc.zipfile
        File? c_mergebam_htmlfile = SE_c_mergebamfqc.htmlfile
        File? c_mergebam_zipfile = SE_c_mergebamfqc.zipfile

        Array[File?]? indv_sp_bam_htmlfile = s_indv_PE_bamfqc.htmlfile
        Array[File?]? indv_sp_bam_zipfile = s_indv_PE_bamfqc.zipfile
        File? sp_mergebam_htmlfile = s_PE_mergebamfqc.htmlfile
        File? sp_mergebam_zipfile = s_PE_mergebamfqc.zipfile

        Array[File?]? indv_cp_bam_htmlfile = c_indv_PE_bamfqc.htmlfile
        Array[File?]? indv_cp_bam_zipfile = c_indv_PE_bamfqc.zipfile
        File? cp_mergebam_htmlfile = c_PE_mergebamfqc.htmlfile
        File? cp_mergebam_zipfile = c_PE_mergebamfqc.zipfile

        File? uno_s_bam_htmlfile = s_uno_PE_bamfqc.htmlfile
        File? uno_s_bam_zipfile = s_uno_PE_bamfqc.zipfile

        File? uno_c_bam_htmlfile = c_uno_PE_bamfqc.htmlfile
        File? uno_c_bam_zipfile = c_uno_PE_bamfqc.zipfile

        #BASICMETRICS
        Array[File?]? s_metrics_out = indv_s_bfs.metrics_out
        Array[File?]? c_metrics_out = indv_c_bfs.metrics_out

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

        Array[File?]? indv_sp_sortedbam = s_indv_PE_mapping.sorted_bam
        Array[File?]? indv_sp_indexbam = s_indv_PE_mapping.bam_index
        Array[File?]? indv_sp_bkbam = s_indv_PE_mapping.bklist_bam
        Array[File?]? indv_sp_bkindexbam = s_indv_PE_mapping.bklist_index
        Array[File?]? indv_sp_rmbam = s_indv_PE_mapping.mkdup_bam
        Array[File?]? indv_sp_rmindexbam = s_indv_PE_mapping.mkdup_index
        Array[File?]? indv_cp_sortedbam = c_indv_PE_mapping.sorted_bam
        Array[File?]? indv_cp_indexbam = c_indv_PE_mapping.bam_index
        Array[File?]? indv_cp_bkbam = c_indv_PE_mapping.bklist_bam
        Array[File?]? indv_cp_bkindexbam = c_indv_PE_mapping.bklist_index
        Array[File?]? indv_cp_rmbam = c_indv_PE_mapping.mkdup_bam
        Array[File?]? indv_cp_rmindexbam = c_indv_PE_mapping.mkdup_index

        File? uno_s_sortedbam = s_uno_PE_mapping.sorted_bam
        File? uno_s_indexstatsbam = s_uno_PE_mapping.bam_index
        File? uno_s_bkbam = s_uno_PE_mapping.bklist_bam
        File? uno_s_bkindexbam = s_uno_PE_mapping.bklist_index
        File? uno_s_rmbam = s_uno_PE_mapping.mkdup_bam
        File? uno_s_rmindexbam = s_uno_PE_mapping.mkdup_index
        File? uno_c_sortedbam = c_uno_PE_mapping.sorted_bam
        File? uno_c_indexstatsbam = c_uno_PE_mapping.bam_index
        File? uno_c_bkbam = c_uno_PE_mapping.bklist_bam
        File? uno_c_bkindexbam = c_uno_PE_mapping.bklist_index
        File? uno_c_rmbam = c_uno_PE_mapping.mkdup_bam
        File? uno_c_rmindexbam = c_uno_PE_mapping.mkdup_index

        File? s_mergebamfile = SE_s_mergebam.mergebam
        File? s_mergebamindex = SE_s_mergeindexstats.indexbam
        File? s_bkbam = SE_s_merge_rmblklist.intersect_out
        File? s_bkindexbam = SE_s_merge_bklist.indexbam
        File? s_rmbam = SE_s_merge_markdup.mkdupbam
        File? s_rmindexbam = SE_s_merge_mkdup.indexbam
        File? c_mergebamfile = SE_c_mergebam.mergebam
        File? c_mergebamindex = SE_c_mergeindexstats.indexbam
        File? c_bkbam = SE_c_merge_rmblklist.intersect_out
        File? c_bkindexbam = SE_c_merge_bklist.indexbam
        File? c_rmbam = SE_c_merge_markdup.mkdupbam
        File? c_rmindexbam = SE_c_merge_mkdup.indexbam

        File? sp_mergebamfile = s_PE_mergebam.mergebam
        File? sp_mergebamindex = s_PE_mergeindexstats.indexbam
        File? sp_bkbam = s_PE_merge_rmblklist.pairtobed_out
        File? sp_bkindexbam = s_PE_merge_bklist.indexbam
        File? sp_rmbam = s_PE_merge_markdup.mkdupbam
        File? sp_rmindexbam = s_PE_merge_mkdup.indexbam
        File? cp_mergebamfile = c_PE_mergebam.mergebam
        File? cp_mergebamindex = c_PE_mergeindexstats.indexbam
        File? cp_bkbam = c_PE_merge_rmblklist.pairtobed_out
        File? cp_bkindexbam = c_PE_merge_bklist.indexbam
        File? cp_rmbam = c_PE_merge_markdup.mkdupbam
        File? cp_rmindexbam = c_PE_merge_mkdup.indexbam

        File? s_fragments_bam = s_fraggraph.fragbamfile
        File? s_fragments_indexbam = s_frag_index.indexbam
        File? c_fragments_bam = c_fraggraph.fragbamfile
        File? c_fragments_indexbam = c_frag_index.indexbam

        #MACS
        File? peakbedfile = SE_macs.peakbedfile
        File? peakxlsfile = SE_macs.peakxlsfile
        File? summitsfile = SE_macs.summitsfile
        File? negativexlsfile = SE_macs.negativepeaks
        File? wigfile = SE_macs.wigfile
        File? ctrlwigfile = SE_macs.ctrlwigfile
        File? all_peakbedfile = SE_all.peakbedfile
        File? all_peakxlsfile = SE_all.peakxlsfile
        File? all_summitsfile = SE_all.summitsfile
        File? all_negativexlsfile = SE_all.negativepeaks
        File? all_wigfile = SE_all.wigfile
        File? all_ctrlwigfile = SE_all.ctrlwigfile
        File? nm_peakbedfile = SE_nomodel.peakbedfile
        File? nm_peakxlsfile = SE_nomodel.peakxlsfile
        File? nm_summitsfile = SE_nomodel.summitsfile
        File? nm_negativexlsfile = SE_nomodel.negativepeaks
        File? nm_wigfile = SE_nomodel.wigfile
        File? nm_ctrlwigfile = SE_nomodel.ctrlwigfile
        File? readme_peaks = SE_addreadme.readme_peaks

        File? only_c_peakbedfile = only_c_macs.peakbedfile
        File? only_c_peakxlsfile = only_c_macs.peakxlsfile
        File? only_c_summitsfile = only_c_macs.summitsfile
        File? only_c_wigfile = only_c_macs.wigfile
        File? only_s_peakbedfile = only_s_macs.peakbedfile
        File? only_s_peakxlsfile = only_s_macs.peakxlsfile
        File? only_s_summitsfile = only_s_macs.summitsfile
        File? only_s_wigfile = only_s_macs.wigfile

        File? only_cp_peakbedfile = only_c_PE_macs.peakbedfile
        File? only_cp_peakxlsfile = only_c_PE_macs.peakxlsfile
        File? only_cp_summitsfile = only_c_PE_macs.summitsfile
        File? only_cp_wigfile = only_c_PE_macs.wigfile
        File? only_sp_peakbedfile = only_s_PE_macs.peakbedfile
        File? only_sp_peakxlsfile = only_s_PE_macs.peakxlsfile
        File? only_sp_summitsfile = only_s_PE_macs.summitsfile
        File? only_sp_wigfile = only_s_PE_macs.wigfile

        File? sp_peakbedfile = PE_macs.peakbedfile
        File? sp_peakxlsfile = PE_macs.peakxlsfile
        File? sp_summitsfile = PE_macs.summitsfile
        File? sp_negativexlsfile = PE_macs.negativepeaks
        File? sp_wigfile = PE_macs.wigfile
        File? sp_ctrlwigfile = PE_macs.ctrlwigfile
        File? sp_all_peakbedfile = PE_all.peakbedfile
        File? sp_all_peakxlsfile = PE_all.peakxlsfile
        File? sp_all_summitsfile = PE_all.summitsfile
        File? sp_all_negativexlsfile = PE_all.negativepeaks
        File? sp_all_wigfile = PE_all.wigfile
        File? sp_all_ctrlwigfile = PE_all.ctrlwigfile
        File? sp_nm_peakbedfile = PE_nomodel.peakbedfile
        File? sp_nm_peakxlsfile = PE_nomodel.peakxlsfile
        File? sp_nm_summitsfile = PE_nomodel.summitsfile
        File? sp_nm_negativexlsfile = PE_nomodel.negativepeaks
        File? sp_nm_wigfile = PE_nomodel.wigfile
        File? sp_nm_ctrlwigfile = PE_nomodel.ctrlwigfile
        File? sp_readme_peaks = PE_addreadme.readme_peaks

        #SICER
        File? scoreisland = SE_sicer.scoreisland
        File? sicer_wigfile = SE_sicer.wigfile
        File? sicer_summary = SE_sicer.summary
        File? sicer_fdrisland = SE_sicer.fdrisland
        File? sp_scoreisland = PE_sicer.scoreisland
        File? sp_sicer_wigfile = PE_sicer.wigfile
        File? sp_sicer_summary = PE_sicer.summary
        File? sp_sicer_fdrisland = PE_sicer.fdrisland

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
        File? supergenes = SE_rose.super_genes
        File? allgenes = SE_rose.all_genes
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
        File? sp_supergenes = PE_rose.super_genes
        File? sp_allgenes = PE_rose.all_genes

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
        File? c_matrices = SE_bamtogff.c_matrices
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
        File? sp_c_matrices = PE_bamtogff.c_matrices
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
        File? c_bigwig = SE_visualization.bigwig
        File? c_norm_wig = SE_visualization.norm_wig
        File? c_tdffile = SE_visualization.tdffile
        File? c_n_bigwig = SE_viznomodel.bigwig
        File? c_n_norm_wig = SE_viznomodel.norm_wig
        File? c_n_tdffile = SE_viznomodel.tdffile
        File? c_a_bigwig = SE_vizall.bigwig
        File? c_a_norm_wig = SE_vizall.norm_wig
        File? c_a_tdffile = SE_vizall.tdffile
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
        File? cp_bigwig = SE_visualization.bigwig
        File? cp_norm_wig = SE_visualization.norm_wig
        File? cp_tdffile = SE_visualization.tdffile
        File? cp_n_bigwig = SE_viznomodel.bigwig
        File? cp_n_norm_wig = SE_viznomodel.norm_wig
        File? cp_n_tdffile = SE_viznomodel.tdffile
        File? cp_a_bigwig = SE_vizall.bigwig
        File? cp_a_norm_wig = SE_vizall.norm_wig
        File? cp_a_tdffile = SE_vizall.tdffile
        File? sp_s_bigwig = PE_vizsicer.bigwig
        File? sp_s_norm_wig = PE_vizsicer.norm_wig
        File? sp_s_tdffile = PE_vizsicer.tdffile

        File? sf_bigwig = s_fraggraph.bigwigfile
        File? sf_tdffile = s_fraggraph.tdffile
        File? sf_wigfile = s_fraggraph.wigfile
        File? cf_bigwig = c_fraggraph.bigwigfile
        File? cf_tdffile = c_fraggraph.tdffile
        File? cf_wigfile = c_fraggraph.wigfile

        #QC-STATS
        Array[File?]? s_qc_statsfile = indv_s_summarystats.statsfile
        Array[File?]? s_qc_htmlfile = indv_s_summarystats.htmlfile
        Array[File?]? s_qc_textfile = indv_s_summarystats.textfile
        File? s_statsfile = SE_s_summarystats.statsfile
        File? s_htmlfile = SE_s_summarystats.htmlfile
        File? s_textfile = SE_s_summarystats.textfile
        Array[File?]? c_qc_statsfile = indv_c_summarystats.statsfile
        Array[File?]? c_qc_htmlfile = indv_c_summarystats.htmlfile
        Array[File?]? c_qc_textfile = indv_c_summarystats.textfile
        File? c_statsfile = SE_c_summarystats.statsfile
        File? c_htmlfile = SE_c_summarystats.htmlfile
        File? c_textfile = SE_c_summarystats.textfile
        File? statsfile = SE_concatstats.statsfile
        File? htmlfile = SE_concatstats.htmlfile
        File? textfile = SE_concatstats.textfile
        File? summaryhtml = SE_overallsummary.summaryhtml
        File? summarytxt = SE_overallsummary.summarytxt

        Array[File?]? sp_qc_statsfile = s_indv_PE_summarystats.statsfile
        Array[File?]? sp_qc_htmlfile = s_indv_PE_summarystats.htmlfile
        Array[File?]? sp_qc_textfile = s_indv_PE_summarystats.textfile
        Array[File?]? cp_qc_statsfile = c_indv_PE_summarystats.statsfile
        Array[File?]? cp_qc_htmlfile = c_indv_PE_summarystats.htmlfile
        Array[File?]? cp_qc_textfile = c_indv_PE_summarystats.textfile
        File? s_qc_mergehtml = s_mergehtmlfile
        File? c_qc_mergehtml = c_mergehtmlfile
        File? s_uno_statsfile = PE_s_summarystats.statsfile
        File? s_uno_htmlfile = PE_s_summarystats.htmlfile
        File? s_uno_textfile = PE_s_summarystats.textfile
        File? c_uno_statsfile = PE_c_summarystats.statsfile
        File? c_uno_htmlfile = PE_c_summarystats.htmlfile
        File? c_uno_textfile = PE_c_summarystats.textfile
        File? sp_statsfile = PE_concatstats.statsfile
        File? sp_htmlfile = PE_concatstats.htmlfile
        File? sp_textfile = PE_concatstats.textfile
        File? s_summaryhtml = PE_merge_s_summarystats.htmlfile
        File? s_summarystats = PE_merge_s_summarystats.statsfile
        File? s_summarytxt = PE_merge_s_summarystats.textfile
        File? c_summaryhtml = PE_merge_c_summarystats.htmlfile
        File? c_summarystats = PE_merge_c_summarystats.statsfile
        File? c_summarytxt = PE_merge_c_summarystats.textfile
        File? sp_summaryhtml = PE_overallsummary.summaryhtml
        File? sp_summarytxt = PE_overallsummary.summarytxt
        File? s_fragsize = s_fraggraph.fragsizepng
        File? c_fragsize = c_fraggraph.fragsizepng
    }
}

