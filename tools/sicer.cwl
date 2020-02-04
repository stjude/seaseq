#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [ sicer ]
class: CommandLineTool

label: SICER new version - Broad Peaks
doc: |
  sicer -t 20190628_KOPTK1-DMSO-MYBL2_AD7124_S12.sorted.bklist.rmdup.bam2bed.bed -s hg19 -egf 0.86 -g 200 -e 100

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing: [ $(inputs.treatmentbedfile) ]

inputs:
  treatmentbedfile:
    type: File
    inputBinding:
      prefix: -t
      position: 1
      valueFrom: $(self.basename)

  controlfile:
    type: File?
    inputBinding:
      prefix: -c
      position: 1

  species: 
    type: string?
    default: "hg19"
    inputBinding:
      prefix: -s
      position: 1
  
  redundancy:
    type: int?
    default: 1
    inputBinding:
      prefix: -rt
      position: 1 

  window:
    type: int?
    default: 200
    inputBinding:
      prefix: -w
      position: 1

  fragment_size:
    type: int?
    default: 150
    inputBinding:
      prefix: -f
      position: 1

  genome_fraction:
    type: double?
    default: 0.86
    inputBinding:
      prefix: -egf
      position: 1

  gapsize:
    type: int?
    default: 200
    inputBinding:
      prefix: -g
      position: 1

  evalue:
    type: int?
    default: 100
    inputBinding:
      prefix: -e
      position: 1

  gzip_wig:
    type: boolean?
    inputBinding:
      position: 998
      shellQuote: false
      prefix: '&& gzip *wig'
    default: true

  outputfolder:
    type: string?
    inputBinding:
      position: 999
      shellQuote: false
      separate: false
      prefix: ' && nameoffolder="'
    default: "SICER_out"

  verifymove:
    type: boolean?
    inputBinding:
      position: 1000 
      shellQuote: false
      prefix: '" && mkdir -p $nameoffolder && mv *W200* $nameoffolder' 
    default: true

outputs:
  sicerDir:
    type: Directory
    outputBinding:
      glob: |
        ${
          return inputs.outputfolder;
        }


