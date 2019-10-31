#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: meme-chip
class: CommandLineTool

label: MEME-ChIP performs comprehensive motif analysis (including motif discovery) 
doc: |
  meme-chip <convert fasta> 

inputs:
  convertfasta:
    type: File
    inputBinding:
      position: 1

outputs:
  outfile:
    type: File
    outputBinding:
      glob: '*/meme.txt'

  htmlfile:
    type: File
    outputBinding:
      glob: '*/meme.html'
