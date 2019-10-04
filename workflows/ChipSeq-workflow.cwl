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


s:name: "ChipSeq-workflow"
s:author:
  - class: s:Person
    s:name: Modupeore Adetunji
    s:email: mailto:madetunj@stjude.org
    s:worksFor:
      - class: s:Organization
        s:name: St.Jude Children's Research Hospital
        s:location: 262 Danny Thomas Place Memphis, TN 38105
        s:department:
          - class: s:Organization
            s:name: Computational Biology
 
s:citation: if any
s:codeRepository: if any
s:dateCreated: "2019-08-01"
s:license: if any
 
doc: |
  Workflow for ChiP-Seq peak calling for single-end fastQ files
  It performs the following steps:
    1. Reference mapping in SAM format using `bowtie`
    2. Convert to BAM, sort, remove duplicates & index are performed using `samtools`.
    3. Peak calling using `macs` & normalization.
    4. Convert normalized peak calls to binary file for visualization using `wigToBigWig`.

  Usage:
    cwlexec -p ChipSeq-workflow.cwl inputparameters.cwl -c LSFconfig.json

$namespaces:
  s: http://schema.org/
 
 
$schemas:
 - https://schema.org/docs/schema_org_rdfa.html
