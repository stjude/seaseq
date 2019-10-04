#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [fastqc]
class: CommandLineTool
label: QC on reads or bam file

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.infile) ]

inputs:
  infile:
    type: File
    inputBinding:
      position: 1
      valueFrom: $(self.basename)
  
outputs:
  htmlfile:
    type: File
    outputBinding: 
      glob: '*fastqc.html'

  zipfile:
    type: File
    outputBinding:
      glob: '*fastqc.zip'
