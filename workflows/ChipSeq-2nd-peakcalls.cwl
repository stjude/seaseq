#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement

inputs:
  reference: Directory
  chromsizes: File
  gfffile: File
  nomodel: boolean?
  wiggle: boolean?
  single_profile: boolean?
  shiftsize: int?
  space: int?
  pvalue: string?
  keep_dup: string?
  motifdatabases: string[]
  flank: int?
  bamfile: File
  zipfile: File
  STATbamoutfile: File
  STATrmdupoutfile: File
  STATbkoutfile: File

outputs:
  peaksbed:
    outputSource: Peaks-Auto/peaksbedfile
    type: File
  
  summits:
    outputSource: Peaks-Auto/summitsfile
    type: File
  
  peaksxls:
    outputSource: Peaks-Auto/peaksxlsfile
    type: File
  
  treatwig:
    outputSource: Peaks-Auto/wigfile
    type: File
  
  statout:
    outputSource: Peaks-Auto/statsfile
    type: File
  
  outBW:
    outputSource: Peaks-Auto/bigwig
    type: File
  
  outtdf:
    outputSource: Peaks-Auto/tdffile
    type: File
  
  ameout:
    type: File
    outputSource: Peaks-Auto/ameout
  
  amehtml:
    type: File
    outputSource: Peaks-Auto/amehtml
  
  memeout:
    type: File
    outputSource: Peaks-Auto/memeout
  
  memehtml:
    type: File
    outputSource: Peaks-Auto/memehtml

  allpeaksbed:
    outputSource: Peaks-All/peaksbedfile
    type: File
  
  allsummits:
    outputSource: Peaks-All/summitsfile
    type: File
  
  allpeaksxls:
    outputSource: Peaks-All/peaksxlsfile
    type: File
  
  alltreatwig:
    outputSource: Peaks-All/wigfile
    type: File
  
  allstatout:
    outputSource: Peaks-All/statsfile
    type: File
  
  alloutBW:
    outputSource: Peaks-All/bigwig
    type: File

  allouttdf:
    outputSource: Peaks-All/tdffile
    type: File

  allameout:
    type: File
    outputSource: Peaks-All/ameout

  allamehtml:
    type: File
    outputSource: Peaks-All/amehtml

  allmemeout:
    type: File
    outputSource: Peaks-All/memeout

  allmemehtml:
    type: File
    outputSource: Peaks-All/memehtml

  sicerbed:
    type: File
    outputSource: SICER/islandbed

  promoterspdf:
    type: File
    outputSource: MetaGene/promoters

  genebodypdf:
    type: File
    outputSource: MetaGene/genebody

  promotersheatmappng:
    type: File
    outputSource: MetaGene/promotersheatmap

  genebodyheatmappng:
    type: File
    outputSource: MetaGene/genebodyheatmap

steps:
  B2Bed:
    in:
      infile: bamfile
    out: [outfile]
    run: ../tools/bamtobed.cwl

  SICER:
    in:
      peaksbedfile: B2Bed/outfile
    out: [islandbed]
    run: ../tools/sicerRB.cwl

  MetaGene:
    in:
      bamfile: bamfile
      gfffile: gfffile
    out: [promoters, genebody, promotersheatmap, genebodyheatmap]
    run: ../tools/bamliquidator.cwl

  Peaks-Auto:
    run: ../subworkflows/peakcalls.cwl
    in:
      bamfile: bamfile
      nomodel: nomodel
      shiftsize: shiftsize
      space: space
      pvalue: pvalue
      wiggle: wiggle
      single_profile: single_profile
      flank: flank
      zipfile: zipfile
      STATbamout: STATbamoutfile
      STATrmdupout: STATrmdupoutfile
      STATbkout: STATbkoutfile
      chromsizes: chromsizes
      reference: reference
      motifdatabases: motifdatabases
    out: [tdffile, peaksbedfile, peaksxlsfile, summitsfile, wigfile, statsfile, bigwig, ameout, amehtml, memeout, memehtml]

  Peaks-All:
    run: ../subworkflows/peakcalls.cwl
    in:
      keep_dup: keep_dup
      bamfile: bamfile
      nomodel: nomodel
      shiftsize: shiftsize
      space: space
      pvalue: pvalue
      wiggle: wiggle
      single_profile: single_profile
      flank: flank
      zipfile: zipfile
      STATbamout: STATbamoutfile
      STATrmdupout: STATrmdupoutfile
      STATbkout: STATbkoutfile
      chromsizes: chromsizes
      reference: reference
      motifdatabases: motifdatabases
    out: [tdffile, peaksbedfile, peaksxlsfile, summitsfile, wigfile, statsfile, bigwig, ameout, amehtml, memeout, memehtml]
