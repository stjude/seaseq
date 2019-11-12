#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement

inputs: 
  nomodel: boolean?
  shiftsize: int?
  space: int?
  pvalue: string?
  keep_dup: string?
  flank: int?
  wiggle: boolean?
  single_profile: boolean?
  bamfile: File
  zipfile: File
  STATbamout: File
  STATrmdupout: File
  STATbkout: File 
  chromsizes: File
  reference: Directory
  motifdatabases: string[]

outputs:
  peaksbedfile:
    type: File
    outputSource: MACS/peaksbedfile

  peaksxlsfile:
    type: File
    outputSource: MACS/peaksxlsfile

  summitsfile:
    type: File
    outputSource: MACS/summitsfile

  wigfile:
    type: File
    outputSource: MACS/wigfile

  statsfile:
    type: File
    outputSource: PeaksQC/statsfile

  bigwig:
    type: File
    outputSource: WigToBig/outfile

  tdffile:
    type: File
    outputSource: toTDF/outfile

  ameout:
    type: File
    outputSource: MOTIFS/ameout

  amehtml:
    type: File
    outputSource: MOTIFS/amehtml

  memeout:
    type: File
    outputSource: MOTIFS/memeout

  memehtml:
    type: File
    outputSource: MOTIFS/memehtml

  summitameout:
    type: File
    outputSource: Summit-MOTIFS/ameout

  summitamehtml:
    type: File
    outputSource: Summit-MOTIFS/amehtml

  summitmemeout:
    type: File
    outputSource: Summit-MOTIFS/memeout

  summitmemehtml:
    type: File
    outputSource: Summit-MOTIFS/memehtml

  rpmwig:
    type: File
    outputSource: RPM/RPMwig

steps:
  MACS:
    in:
      treatmentfile: bamfile
      nomodel: nomodel
      shiftsize: shiftsize
      space: space
      pvalue: pvalue
      keep_dup: keep_dup
      wiggle: wiggle
      single_profile: single_profile
    out: [peaksbedfile, peaksxlsfile, summitsfile, wigfile]
    run: ../tools/macs1call.cwl

  PeaksQC:
    in:
      fastqczip: zipfile
      bamfile: bamfile
      peaksbed: MACS/peaksbedfile
      bamflag: STATbamout
      rmdupflag: STATrmdupout
      bkflag: STATbkout
    out: [ statsfile ]
    run: ../tools/summarystats.cwl

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

  toTDF:
    in:
      wigfile: MACS/wigfile
    out: [outfile]
    run: ../tools/igvtdf.cwl

  FlankBED:
    in:
      bedfile: MACS/summitsfile
      flank: flank
    out: [outfile]
    run: ../tools/flankbed.cwl

  MOTIFS:
    run: ../subworkflows/motifs.cwl
    in:
      reference: reference
      bedfile: MACS/peaksbedfile
      motifdatabases: motifdatabases
    out: [ameout, amehtml, memeout, memehtml]

  Summit-MOTIFS:
    run: ../subworkflows/motifs.cwl
    in:
      reference: reference
      bedfile: FlankBED/outfile
      motifdatabases: motifdatabases
    out: [ameout, amehtml, memeout, memehtml]


doc: |
  Workflow performs Peak Calls from `bamfile` input using MACSv1 and SICERv1.1 on single-end (SE) reads.
  In addition MOTIFs (Enriched and Discovered) are called using the MEME-suite (AME & MEMe-chip), and 
  METAGENES are also called using bamliquidator, i.e. bamtoGFF.py v1.
