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
  redundancy: int?
  window: int?
  fragment_size: int?
  genome_fraction: double?
  gapsize: int?
  evalue: int?

  # ROSE
  feature: string?

# addendum from ChipSeq-1st-mapping step
  rmdupbamfile: 
    type: File
    secondaryFiles:
      - .bai
  bklistbamfile: 
    type: File
    secondaryFiles:
      - .bai
  readzipfile: File
  STATbamoutfile: File
  STATrmdupoutfile: File
  STATbkoutfile: File
  fastqmetricsfile: File

outputs:
# MACS-AUTO
  macsDir:
    type: Directory
    outputSource: MACS-Auto/macsDir

# MACS-ALL
  allmacsDir:
    type: Directory
    outputSource: MACS-All/macsDir

# MACS-NM
  nmmacsDir:
    type: Directory
    outputSource: MACS-NM/macsDir

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

# MOTIFs & Summits output
  bedfasta:
    type: File
    outputSource: MOTIFS/bedfasta

  flankbed:
    type: File
    outputSource: FlankBED/outfile
    
  memechipdir:
    type: Directory
    outputSource: MOTIFS/memechipdir

  summitmemechipdir:
    type: Directory
    outputSource: SummitMOTIFS/memechipdir

  amedir:
    type: Directory
    outputSource: MOTIFS/amedir

  summitamedir:
    type: Directory
    outputSource: SummitMOTIFS/amedir
    
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

  textfile:
    type: File
    outputSource: PeaksQC/textfile

steps:
# PEAK CALLING & VISUALS
  MACS-Auto:
    requirements:
      ResourceRequirement:
        ramMax: 10000
        coresMin: 1
    in:
      treatmentfile: bklistbamfile
      space: space
      pvalue: pvalue
      wiggle: wiggle
      single_profile: single_profile
    out: [ peaksbedfile, peaksxlsfile, summitsfile, wigfile, macsDir ]
    run: ../tools/macs1call.cwl

  WIG-Auto:
    in:
      wigfile: MACS-Auto/wigfile
      peaksxls: MACS-Auto/peaksxlsfile
      chromsizes: chromsizes
    out: [ rpmwig, outBW, outtdf ]
    run: ../subworkflows/visualization.cwl

  MACS-All:
    requirements:
      ResourceRequirement:
        ramMax: 10000
        coresMin: 1
    in:
      treatmentfile: bklistbamfile
      keep_dup: keep_dup
      space: space
      pvalue: pvalue
      wiggle: wiggle
      single_profile: single_profile
    out: [ peaksbedfile, peaksxlsfile, summitsfile, wigfile, macsDir ]
    run: ../tools/macs1call.cwl

  WIG-All:
    in:
      wigfile: MACS-All/wigfile
      peaksxls: MACS-All/peaksxlsfile
      chromsizes: chromsizes
    out: [ rpmwig, outBW, outtdf ]
    run: ../subworkflows/visualization.cwl

  MACS-NM:
    requirements:
      ResourceRequirement:
        ramMax: 10000
        coresMin: 1
    in:
      treatmentfile: bklistbamfile
      space: space
      wiggle: wiggle
      single_profile: single_profile
    out: [ peaksbedfile, peaksxlsfile, summitsfile, wigfile, macsDir ]
    run: ../tools/macs1nm.cwl

  WIG-NM:
    in:
      wigfile: MACS-NM/wigfile
      peaksxls: MACS-NM/peaksxlsfile
      chromsizes: chromsizes
    out: [ rpmwig, outBW, outtdf ]
    run: ../subworkflows/visualization.cwl

# MOTIF analysis
  MOTIFS:
    in:
      reference: reference
      bedfile: MACS-Auto/peaksbedfile
      motifdatabases: motifdatabases
    out: [memechipdir, amedir, bedfasta]
    run: ../subworkflows/motifs.cwl

#SUMMIT-MOTIF analysis
  FlankBED:
    in:
      bedfile: MACS-Auto/summitsfile
      flank: flank
    out: [outfile]
    run: ../tools/flankbed.cwl
  
  SummitMOTIFS:
    in:
      reference: reference
      bedfile: FlankBED/outfile
      motifdatabases: motifdatabases
    out: [memechipdir, amedir, bedfasta]
    run: ../subworkflows/motifs.cwl
    
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
    requirements:
      ResourceRequirement:
        ramMax: 10000
        coresMin: 1
    in:
      species: species
      redundancy: redundancy
      window: window
      fragment_size: fragment_size
      genome_fraction: genome_fraction
      gapsize: gapsize
      evalue: evalue
      treatmentbedfile: B2Bed/outfile
    run: ../tools/sicer.cwl
    out: [ sicerDir ]

# ROSE enhancer caller
  ROSE:
    requirements:
      ResourceRequirement:
        ramMax: 20000
        coresMin: 1
    in:
      species: species
      feature: feature
      gtffile: gtffile
      bamfile: bklistbamfile
      fileA: MACS-All/peaksbedfile
      fileB: MACS-Auto/peaksbedfile
    out: [ RoseDir ]
    run: ../tools/roseNC.cwl

# Quality Control & Statistics
  Bklist2Bed:
    in:
      infile: bklistbamfile
    out: [ outfile ]
    run: ../tools/bamtobed.cwl

  SortBed:
    requirements:
      ResourceRequirement:
        ramMax: 10000
        coresMin: 1
    in:
      infile: Bklist2Bed/outfile
    out: [outfile]
    run: ../tools/sortbed.cwl

  runSPP:
    requirements:
      ResourceRequirement:
        ramMax: 10000
        coresMin: 1
    in:
      infile: bklistbamfile
    out: [spp_out]
    run: ../tools/runSPP.cwl

  CountIntersectBed:
    requirements:
      ResourceRequirement:
        ramMax: 10000
        coresMin: 1
    in:
      peaksbed: MACS-Auto/peaksbedfile
      bamtobed: SortBed/outfile
    out: [outfile]
    run: ../tools/intersectbed.cwl

  PeaksQC:
    requirements:
      ResourceRequirement:
        ramMax: 10000
        coresMin: 1
    in:
      fastqmetrics: fastqmetricsfile
      fastqczip: readzipfile
      sppfile: runSPP/spp_out
      bambed: Bklist2Bed/outfile
      countsfile: CountIntersectBed/outfile
      peaksxls: MACS-Auto/peaksxlsfile
      bamflag: STATbamoutfile
      rmdupflag: STATrmdupoutfile
      bkflag: STATbkoutfile
      rosedir: ROSE/RoseDir
    out: [ statsfile, htmlfile, textfile ]
    run: ../tools/summarystatsv2.cwl
