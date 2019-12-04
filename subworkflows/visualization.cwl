#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement

inputs: 
  chromsizes: File
  wigfile: File
  peaksxls: File
  

outputs:
  rpmwig:
    type: File
    outputSource: RPM/RPMwig

  outBW:
    outputSource: WIG/outfile
    type: File

  outtdf:
    outputSource: TDF/outfile
    type: File


steps:
  RPM:
    in:
      wigfile: wigfile
      peaksxls: peaksxls
    out: [RPMwig]
    run: ../tools/normalize.cwl

  WIG:
    in:
      infile: RPM/RPMwig
      chromsizes: chromsizes
    out: [outfile]
    run: ../tools/wigtobigwig.cwl

  TDF:
    in:
      wigfile: wigfile
    out: [outfile]
    run: ../tools/igvtdf.cwl

