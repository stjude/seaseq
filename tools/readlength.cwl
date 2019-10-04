#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [readLength.sh]
class: CommandLineTool

inputs:
  fastqfile:
    type: File
    label: "FastQfiles"
    inputBinding:
      position: 1

outputs:
  readLength:
    type: File
    label: random file containing read length
    outputBinding:
      glob: '*.rdl'
