#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [ sicer ]
class: CommandLineTool

label: SICER new version - Broad Peaks
doc: |
  sicer -t 20190628_KOPTK1-DMSO-MYBL2_AD7124_S12.sorted.bklist.rmdup.bam2bed.bed -s hg19 -egf 0.86 -g 200 -e 100

requirements:
- class: InitialWorkDirRequirement
  listing: [ $(inputs.treatmentbedfile) ]

inputs:
  treatmentbedfile:
    type: File
    inputBinding:
      prefix: -t
      valueFrom: $(self.basename)

  controlfile:
    type: File?
    inputBinding:
      prefix: -c

  species: 
    type: string?
    default: "hg19"
    inputBinding:
      prefix: -s
  
  redundancy:
    type: int?
    default: 1
    inputBinding:
      prefix: -rt 

  window:
    type: int?
    default: 200
    inputBinding:
      prefix: -w

  fragment_size:
    type: int?
    default: 150
    inputBinding:
      prefix: -f

  genome_fraction:
    type: double?
    default: 0.86
    inputBinding:
      prefix: -egf

  gapsize:
    type: int?
    default: 200
    inputBinding:
      prefix: -g

  evalue:
    type: int?
    default: 100
    inputBinding:
      prefix: -e

  outfile_txt:
    type: string?
    inputBinding:
      position: 999
      shellQuote: false
      prefix: '&& gzip *wig'
    default: ""

  outdir:
    type: string?
    inputBinding:
      position: 1000
      shellQuote: false
      prefix: '&& mkdir -p SICER_out && mv *W200* SICER_out'
    default: ""

outputs:
  sicerDir:
    type: Directory
    outputBinding:
      glob: "SICER_out"
