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
    label: "GTF file"
    inputBinding:
      position: 1

  bamfile:
    type: File
    label: "BAM file"
    inputBinding:
      position: 2
#    secondaryFiles: 
#      - .bai
    
  outputdir:
    type: string
    label: "Output directory name"
    default: "ROSE_out"
    inputBinding:
      position: 3

  feature:
    type: string?
    label: "Feature Type"
    default: "gene"
    inputBinding:
      position: 4

  species:
    type: string?
    label: "Genome name"
    default: "hg19"
    inputBinding:
      position: 5

  fileA:
    type: File
    label: "MACS Auto BED file"
    inputBinding:
      position: 6

  fileB:
    type: File
    label: "MACS All BED file"
    inputBinding:
      position: 7

outputs:
  RoseDir:
    type: Directory
    label: "ROSE output directory"
    outputBinding:
      glob: |
        ${
          return inputs.outputdir;
        }
