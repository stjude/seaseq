#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ normalize_WIG_to_RPM.pl ]
doc: |
  Normalize peaks wig score per million


hints:
  DockerRequirement:
    dockerPull: madetunj/seaseq:v0.0.1


inputs:
  wigfile:
    type: File
    label: "Wig file"
    inputBinding:
      position: 1
  
  peaksxls:
    type: File
    label: "MACS Peaks xls file"
    inputBinding:
      position: 2

    
outputs:
  RPMwig:
    type: File
    label: "normalized wig scores"
    outputBinding:
      glob: '*.RPM.wig.gz'
