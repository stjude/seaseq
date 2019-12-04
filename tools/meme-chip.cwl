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
  outDir:
    type: Directory
    outputBinding:
      glob: "memechip_out"

  memeDir:
    type: Directory
    outputBinding:
      glob: "meme_out"
