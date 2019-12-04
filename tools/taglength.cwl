#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [tagLength.sh]
class: CommandLineTool

inputs:
  datafile:
    type: File
    inputBinding:
      position: 1

outputs:
  tagLength:
    type: File
    label: random file containing mean read length
    outputBinding:
      glob: '*.rdl'
