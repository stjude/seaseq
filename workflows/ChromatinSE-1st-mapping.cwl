#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement

inputs:
  reference: Directory
  fastqfile: File
  chromsizes: File
  blacklistfile: File
  best_alignments: boolean?
  good_alignments: int?
  limit_alignments: int?
  processors: int?

outputs:
  sam_sort:
    outputSource: SamSort/outfile
    type: File

  fastq_metrics:
    outputSource: BasicMetrics/metrics_out
    type: File

  rmdup_bam:
    outputSource: SamRMDup/outfile
    type: File

  rmdup_index:
    outputSource: SamIndex/outfile
    type: File

  bklist_bam:
    outputSource: BkList/outfile
    type: File

  bklist_index: 
    outputSource: BkIndex/outfile
    type: File

  bamqc_html:
    outputSource: BamQC/htmlfile
    type: File

  bamqc_zip:
    outputSource: BamQC/zipfile
    type: File

  readqc_zip:
    outputSource: ReadQC/zipfile
    type: File

  readqc_html:
    outputSource: ReadQC/htmlfile
    type: File

  stat_bk:
    outputSource: STATbk/outfile
    type: File

  stat_bam:
    outputSource: STATbam/outfile
    type: File

  stat_rmdup:
    outputSource: STATrmdup/outfile
    type: File
    

steps:
  BasicMetrics:
    in: 
      fastqfile: fastqfile
    out: [metrics_out]
    run: ../tools/basicfastqstats.cwl

  TagLen:
    in: 
      datafile: BasicMetrics/metrics_out
    out: [tagLength]
    run: ../tools/taglength.cwl
   
  ReadQC:
    in:
      infile: fastqfile
    out: [htmlfile, zipfile]
    run: ../tools/fastqc.cwl

  Bowtie:
    run: ../tools/bowtie.cwl
    in:
      readLengthFile: TagLen/tagLength
      best_alignments: best_alignments
      good_alignments: good_alignments
      fastqfile: fastqfile
      limit_alignments: limit_alignments
      processors: processors
      reference: reference
    out: [samfile]

  SamView:
    in:
      infile: Bowtie/samfile
    out: [outfile]
    run: ../tools/samtools-view.cwl

  BamQC:
    in:
      infile: SamView/outfile
    out: [htmlfile, zipfile]
    run: ../tools/fastqc.cwl

  SamSort:
    in:
      infile: SamView/outfile
    out: [outfile]
    run: ../tools/samtools-sort.cwl

  BkList:
    in:
      infile: SamSort/outfile
      blacklistfile: blacklistfile
    out: [outfile]
    run: ../tools/blacklist.cwl

  BkIndex:
    in:
      infile: BkList/outfile
    out: [bam2file, outfile]
    run: ../tools/samtools-index.cwl

  SamRMDup:
    in:
      infile: BkList/outfile
    out: [outfile]
    run: ../tools/samtools-mkdupr.cwl

  SamIndex:
    in:
      infile: SamRMDup/outfile
    out: [bam2file, outfile]
    run: ../tools/samtools-index.cwl

  STATbam:
    in:
      infile: SamView/outfile
    out: [outfile]
    run: ../tools/samtools-flagstat.cwl

  STATrmdup:
    in:
      infile: SamRMDup/outfile
    out: [outfile]
    run: ../tools/samtools-flagstat.cwl

  STATbk:
    in:
      infile: BkList/outfile
    out: [outfile]
    run: ../tools/samtools-flagstat.cwl


doc: |
  Runs ChIP-Seq SE Mapping FastQ SE files to generate BAM file for step 2 in ChIP-Seq Pipeline.
