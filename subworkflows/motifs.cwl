#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement

inputs: 
  reference: Directory
  bedfile: File
  motifdatabases: string[]

outputs:
  memechipoutdir:
    type: Directory
    outputSource: MEMECHIP/outDir

  memeoutdir:
    type: Directory
    outputSource: MEMECHIP/memeDir

  ameoutdir:
    type: Directory
    outputSource: AME/outDir

steps:
  MEMECHIP:
    run: ../tools/meme-chip.cwl
    in:
      convertfasta: BEDfasta/outfile
    out: [outDir, memeDir]

  AME:
    run: ../tools/ame.cwl
    in:
      convertfasta: BEDfasta/outfile
      motifdatabases: motifdatabases
    out: [outDir]

  BEDfasta:
    in:
      reference: reference
      bedfile: bedfile
    out: [outfile]
    run: ../tools/bedfasta.cwl


doc: |
  Workflow calls MOTIFS both enriched and discovered using the MEME-suite (AME & MEME-chip).
