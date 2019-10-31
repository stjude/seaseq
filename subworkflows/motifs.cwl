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
  ameout:
    type: File
    outputSource: AME/outfile

  amehtml:
    type: File
    outputSource: AME/htmlfile

  memeout:
    type: File
    outputSource: MEME/outfile

  memehtml:
    type: File
    outputSource: MEME/htmlfile

steps:
  BEDfasta:
    in:
      reference: reference
      bedfile: bedfile
    out: [outfile]
    run: ../tools/bedfasta.cwl

  AME:
    in:
      convertfasta: BEDfasta/outfile
      motifdatabases: motifdatabases
    out: [outfile, htmlfile]
    run: ../tools/ame.cwl

  MEME:
    in:
      convertfasta: BEDfasta/outfile
    out: [outfile, htmlfile]
    run: ../tools/meme-chip.cwl


doc: |
  Workflow calls MOTIFS both enriched and discovered using the MEME-suite (AME & MEME-chip).
