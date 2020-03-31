#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [ROSE-local.sh]
class: CommandLineTool
label: ROSE - calling Enhancers and Super-enhancers
doc: |
  ROSE_call.sh <gtf file> <bam file> ROSE_out genes hg19 <bed file 1> <bed file 2>

requirements:
- class: InlineJavascriptRequirement

inputs:
  gtffile:
    type: File
    inputBinding:
      position: 1

  bamfile:
    type: File
    inputBinding:
      position: 2
#    secondaryFiles: 
#      - .bai
    
  outputdir:
    type: string
    default: "ROSE_out"
    inputBinding:
      position: 3

  feature:
    type: string?
    default: "gene"
    inputBinding:
      position: 4

  species:
    type: string?
    default: "hg19"
    inputBinding:
      position: 5

  fileA:
    type: File
    inputBinding:
      position: 6

  fileB:
    type: File
    inputBinding:
      position: 7

outputs:
  RoseDir:
    type: Directory
    outputBinding:
      glob: |
        ${
          return inputs.outputdir;
        }
