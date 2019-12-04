#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: ame
class: CommandLineTool

label: AME - Analysis of Motif Enrichment
doc: |
  ame <convert fasta> <motif-databases>


inputs:
  convertfasta:
    type: File
    inputBinding:
      position: 1
  
  motifdatabases:
    type: 
      type: array
      items: string
    inputBinding:
      position: 2


outputs:
  outDir:
    type: Directory
    outputBinding:
      glob: "ame_out"

