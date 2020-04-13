#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

doc: |
  Create visualization files: TDF for IGV and BigWig for UCSC

requirements:
  - class: SubworkflowFeatureRequirement

inputs: 
  chromsizes: 
    type: File
    label: "Chromosome sizes tab file"

  wigfile:
    type: File
    label: "Wig file"

  peaksxls: 
    type: File
    label: "MACS1 Peaks xls file"


outputs:
  rpmwig:
    type: File
    outputSource: RPM/RPMwig
    label: "normalized wig scores"

  outBW:
    outputSource: WIG/outfile
    type: File
    label: "BigWig file"

  outtdf:
    outputSource: TDF/outfile
    type: File
    label: "IGV TDF file"


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

