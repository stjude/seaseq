#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand: [ normalize_WIG_to_RPM.pl ]
inputs:
  wigfile:
    type: File
    inputBinding:
      position: 1
  
  peaksxls:
    type: File
    inputBinding:
      position: 2
    
outputs:
  RPMwig:
    type: File
    outputBinding:
      glob: '*.RPM.wig.gz'
