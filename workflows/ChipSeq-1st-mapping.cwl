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
  bam:
    outputSource: SamRMDup/outfile
    type: File

  index:
    outputSource: SamIndex/outfile
    type: File

  bamqchtml:
    outputSource: BamQC/htmlfile
    type: File

  bamqczip:
    outputSource: BamQC/zipfile
    type: File

  readqczip:
    outputSource: ReadQC/zipfile
    type: File

  readqchtml:
    outputSource: ReadQC/htmlfile
    type: File

  statbk:
    outputSource: STATbk/outfile
    type: File

  statbam:
    outputSource: STATbam/outfile
    type: File

  statrmdup:
    outputSource: STATrmdup/outfile
    type: File
    
steps:
  ReadLen:
    in:
      fastqfile: fastqfile
    out: [readLength]
    run: ../tools/readlength.cwl
   

  ReadQC:
    in:
      infile: fastqfile
    out: [htmlfile, zipfile]
    run: ../tools/fastqc.cwl

  Bowtie:
    run: ../tools/bowtie.cwl
    in:
      readLengthFile: ReadLen/readLength
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
      infile: SamRMDup/outfile
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

  SamRMDup:
    in:
      infile: BkList/outfile
    out: [outfile]
    run: ../tools/samtools-rmdup.cwl

  SamIndex:
    in:
      infile: SamRMDup/outfile
    out: [bam2file, outfile]
    run: ../tools/samtools-index.cwl

  STATbam:
    in:
      infile: SamView/outfile
    out: [outfile]
    run: ../tools/sambamba.cwl

  STATrmdup:
    in:
      infile: SamRMDup/outfile
    out: [outfile]
    run: ../tools/sambamba.cwl

  STATbk:
    in:
      infile: BkList/outfile
    out: [outfile]
    run: ../tools/sambamba.cwl

doc: |
  Runs ChIP-Seq SE Mapping FastQ SE files to generate BAM file for step 2 in ChIP-Seq Pipeline.
