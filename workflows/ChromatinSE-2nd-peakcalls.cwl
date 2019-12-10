#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement


inputs:
# main files & directorys
  reference: Directory
  chromsizes: File
  gtffile: File
  motifdatabases: string[]

# optional parameters
  # MACS
  nomodel: boolean?
  wiggle: boolean?
  single_profile: boolean?
  shiftsize: int?
  space: int?
  pvalue: string?
  keep_dup: string?
  flank: int?
  
  # SICER & ROSE
  species: string?
  
  # SICER
  redundancy: string?
  window: int?
  fragment_size: int?
  genome_fraction: double?
  gapsize: int?
  evalue: int?

  # ROSE
  feature: string?

# addendum from ChipSeq-1st-mapping step
  rmdupbamfile: File
  bklistbamfile: File
  bklistbamindex: File
  readzipfile: File
  STATbamoutfile: File
  STATrmdupoutfile: File
  STATbkoutfile: File
  fastqmetricsfile: File

outputs:
# MACS-AUTO
  peaksbed:
    outputSource: MACS-Auto/peaksbedfile
    type: File

  summits:
    outputSource: MACS-Auto/summitsfile
    type: File

  peaksxls:
    outputSource: MACS-Auto/peaksxlsfile
    type: File

  treatwig:
    outputSource: MACS-Auto/wigfile
    type: File

  modelR:
    type: File
    outputSource: MACS-Auto/modelR

# MACS-ALL
  allpeaksbed:
    outputSource: MACS-All/peaksbedfile
    type: File

  allsummits:
    outputSource: MACS-All/summitsfile
    type: File

  allpeaksxls:
    outputSource: MACS-All/peaksxlsfile
    type: File

  alltreatwig:
    outputSource: MACS-All/wigfile
    type: File

  allmodelR:  
    type: File  
    outputSource: MACS-All/modelR  

# MACS-NM
  nmpeaksbed:
    outputSource: MACS-NM/peaksbedfile
    type: File

  nmsummits:
    outputSource: MACS-NM/summitsfile
    type: File

  nmpeaksxls:
    outputSource: MACS-NM/peaksxlsfile
    type: File

  nmtreatwig:
    outputSource: MACS-NM/wigfile
    type: File

#VISUAL
  rpmwig:
    type: File
    outputSource: WIG-Auto/rpmwig

  outBW:
    outputSource: WIG-Auto/outBW
    type: File

  outtdf:
    outputSource: WIG-Auto/outtdf
    type: File

  allrpmwig:
    type: File
    outputSource: WIG-All/rpmwig

  alloutBW:
    outputSource: WIG-All/outBW
    type: File

  allouttdf:
    outputSource: WIG-All/outtdf
    type: File

  nmrpmwig:
    type: File
    outputSource: WIG-NM/rpmwig

  nmoutBW:
    outputSource: WIG-NM/outBW
    type: File

  nmouttdf:
    outputSource: WIG-NM/outtdf
    type: File

# MOTIFs output
  bedfasta:
    type: File
    outputSource: BEDfasta/outfile

  memechipdir:
    type: Directory
    outputSource: MEMECHIP/outDir

  ameoutdir:
    type: Directory
    outputSource: AME/outDir

# METAGENE output
  metagenesDir:
    type: Directory
    outputSource: MetaGene/metagenesDir

# SICER output
  sicerDir:
    type: Directory
    outputSource: SICER/sicerDir

# ROSE output
  roseoutput:
    type: Directory
    outputSource: ROSE/RoseDir

# QC Control & Statistics output
  statsfile:
    type: File
    outputSource: PeaksQC/statsfile

  htmlfile:
    type: File
    outputSource: PeaksQC/htmlfile


steps:
# PEAK CALLING & VISUALS
  MACS-Auto:
    in:
      treatmentfile: bklistbamfile
      space: space
      pvalue: pvalue
      wiggle: wiggle
      single_profile: single_profile
    out: [ peaksbedfile, peaksxlsfile, summitsfile, wigfile, modelR ]
    run: ../tools/macs1call.cwl

  WIG-Auto:
    in:
      wigfile: MACS-Auto/wigfile
      peaksxls: MACS-Auto/peaksxlsfile
      chromsizes: chromsizes
    out: [ rpmwig, outBW, outtdf ]
    run: ../subworkflows/visualization.cwl

  MACS-All:
    in:
      treatmentfile: bklistbamfile
      keep_dup: keep_dup
      space: space
      pvalue: pvalue
      wiggle: wiggle
      single_profile: single_profile
    out: [ peaksbedfile, peaksxlsfile, summitsfile, wigfile, modelR ]
    run: ../tools/macs1call.cwl

  WIG-All:
    in:
      wigfile: MACS-All/wigfile
      peaksxls: MACS-All/peaksxlsfile
      chromsizes: chromsizes
    out: [ rpmwig, outBW, outtdf ]
    run: ../subworkflows/visualization.cwl

  MACS-NM:
    in:
      treatmentfile: bklistbamfile
      space: space
      pvalue: pvalue
      wiggle: wiggle
      single_profile: single_profile
    out: [ peaksbedfile, peaksxlsfile, summitsfile, wigfile ]
    run: ../tools/macs1nm.cwl

  WIG-NM:
    in:
      wigfile: MACS-NM/wigfile
      peaksxls: MACS-NM/peaksxlsfile
      chromsizes: chromsizes
    out: [ rpmwig, outBW, outtdf ]
    run: ../subworkflows/visualization.cwl

# MOTIF analysis
  BEDfasta:
    in:
      reference: reference
      bedfile: MACS-Auto/peaksbedfile
    out: [ outfile ]
    run: ../tools/bedfasta.cwl

  MEMECHIP:
    run: ../tools/meme-chip.cwl
    in:
      convertfasta: BEDfasta/outfile
    out: [ outDir ]

  AME:
    run: ../tools/ame.cwl
    in:
      convertfasta: BEDfasta/outfile
      motifdatabases: motifdatabases
    out: [ outDir ]

# METAGENE analysis
  MetaGene:
    in:
      bamfile: rmdupbamfile
      gtffile: gtffile
      chromsizes: chromsizes
    out:  [ metagenesDir ] 
    run: ../tools/bamtogff.cwl

# SICER broad peaks caller
  B2Bed:
    in:
      infile: rmdupbamfile
    out: [ outfile ]
    run: ../tools/bamtobed.cwl

  SICER:
    in:
      species: species
      redundancy: redundancy
      window: window
      fragment_size: fragment_size
      genome_fraction: genome_fraction
      gapsize: gapsize
      evalue: evalue
      peaksbedfile: B2Bed/outfile
    run: ../tools/sicerRB.cwl
    out: [ sicerDir ]

# ROSE enhancer caller
  ROSE:
    in:
      species: species
      feature: feature
      gtffile: gtffile
      bamfile: bklistbamfile
      bamindex: bklistbamindex
      fileA: MACS-All/peaksbedfile
      fileB: MACS-Auto/peaksbedfile
    out: [ RoseDir ]
    run: ../tools/roseNC.cwl

# Quality Control & Statistics
  PeaksQC:
    in:
      fastqmetrics: fastqmetricsfile
      fastqczip: readzipfile
      bamfile: bklistbamfile
      peaksbed: MACS-Auto/peaksbedfile
      peaksxls: MACS-Auto/peaksxlsfile
      bamflag: STATbamoutfile
      rmdupflag: STATrmdupoutfile
      bkflag: STATbkoutfile
      rosedir: ROSE/RoseDir
    out: [ statsfile, htmlfile ]
    run: ../tools/summarystats.cwl
