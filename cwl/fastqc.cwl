#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [fastqc]
class: CommandLineTool
label: QC on reads or bam file

hints:
  DockerRequirement:
    dockerPull: madetunj/fastqc:v0.11.9

requirements:
- class: InitialWorkDirRequirement
  listing: [ $(inputs.infile) ]

inputs:
  infile:
    type: File
    label: "Input file"
    inputBinding:
      position: 1
  
  outputdirectory:
    type: string?
    label: "Output directory location"
    inputBinding:
      position: 2
      prefix: '-o'
    default: './'

outputs:
  htmlfile:
    type: File
    label: "FastQC HTML file"
    outputBinding: 
      glob: '*_fastqc.html'

  zipfile:
    type: File
    label: "FastQC ZIP file"
    outputBinding:
      glob: '*_fastqc.zip'
