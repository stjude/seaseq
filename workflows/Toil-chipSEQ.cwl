#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement

inputs:
  reference: Directory
  fastqfile: File
  chromsizes: File
  output_prefix: string?
  best_alignments: boolean?
  nomodel: boolean?
  wiggle: boolean?
  single_profile: boolean?
  good_alignments: int?
  limit_alignments: int?
  processors: int?
  shiftsize: int?
  space: int?


outputs:
  sam:
    outputSource: Bowtie/samfile
    type: File

  bam: 
    outputSource: SamView/outfile
    type: File

  sort:
    outputSource: SamSort/outfile
    type: File

  rmdup:
    outputSource: SamRMDup/outfile
    type: File

  index:
    outputSource: SamIndex/outfile
    type: File
  
  peaksbed:
    outputSource: MACS/peaksbedfile
    type: File

  summits:
    outputSource: MACS/summitsfile
    type: File

  peaksxls:
    outputSource: MACS/peaksxlsfile
    type: File

  treat:
    outputSource: MACS/wigfile
    type: File

  outRPM:
    outputSource: RPM/RPMwig
    type: File

  outBW:
    outputSource: WigToBig/outfile
    type: File


steps:
  ReadLen:
    in: 
      fastqfile: fastqfile
    out: [readLength]
    run: ../tools/readlength.cwl

  Bowtie:
    requirements:
      ResourceRequirement:
        coresMin: 20
    run: ../tools/bowtie.cwl
    in:
      output_prefix: output_prefix
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

  SamSort:
    in:
      infile: SamView/outfile
    out: [outfile]
    run: ../tools/samtools-sort.cwl

  SamRMDup:
    in:
      infile: SamSort/outfile
    out: [outfile]
    run: ../tools/samtools-rmdup.cwl

  SamIndex:
    in:
      infile: SamRMDup/outfile
    out: [bam2file, outfile]
    run: ../tools/samtools-index.cwl

  MACS:
    in:
      treatmentfile: SamIndex/bam2file
      nomodel: nomodel
      shiftsize: shiftsize
      space: space
      wiggle: wiggle
      single_profile: single_profile
    out: [peaksbedfile, peaksxlsfile, summitsfile, wigfile]
    run: ../tools/macs.cwl

  RMDupSort:
    in:
      infile: SamRMDup/outfile
    out: [outfile]
    run: ../tools/samtools-sort.cwl

  RMDupIndex:
    in:
      infile: RMDupSort/outfile
    out: [bam2file, outfile]
    run: ../tools/samtools-index.cwl

  RPM:
    in:
      wigfile: MACS/wigfile
      peaksxls: MACS/peaksxlsfile
    out: [RPMwig]
    run: ../tools/normalize.cwl

  WigToBig:
    in:
      infile: RPM/RPMwig
      chromsizes: chromsizes
    out: [outfile]
    run: ../tools/wigtobigwig.cwl
