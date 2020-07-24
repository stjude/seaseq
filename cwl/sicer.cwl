#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [ sicer ]
class: CommandLineTool
label: SICER new version - Broad Peaks
doc: |
  sicer -t 20190628_KOPTK1-DMSO-MYBL2_AD7124_S12.sorted.bklist.rmdup.bam2bed.bed -s hg19 -egf 0.86 -g 200 -e 100


hints:
  DockerRequirement:
    dockerPull: madetunj/sicer:v1.0.2 


requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing: [ $(inputs.treatmentbedfile) ]


inputs:
  treatmentbedfile:
    type: File
    label: "BED file"
    inputBinding:
      prefix: -t
      position: 1
      valueFrom: $(self.basename)

  controlfile:
    type: File?
    label: "Control BED file"
    inputBinding:
      prefix: -c
      position: 1

  species: 
    type: string?
    label: "Genome name"
    default: "hg19"
    inputBinding:
      prefix: -s
      position: 1
  
  redundancy:
    type: int?
    label: "Number of identical reads allowed"
    default: 1
    inputBinding:
      prefix: -rt
      position: 1 

  window:
    type: int?
    label: "Resolution of SICER"
    default: 200
    inputBinding:
      prefix: -w
      position: 1

  fragment_size:
    type: int?
    label: "Fragment shift size"
    default: 150
    inputBinding:
      prefix: -f
      position: 1

  genome_fraction:
    type: double?
    label: "Effectve Genome Fraction"
    default: 0.86
    inputBinding:
      prefix: -egf
      position: 1

  gapsize:
    type: int?
    label: "Gap size"
    default: 200
    inputBinding:
      prefix: -g
      position: 1

  evalue:
    type: int?
    label: "E-value"
    default: 100
    inputBinding:
      prefix: -e
      position: 1

  gzip_wig:
    type: boolean?
    label: "GZIP wig files"
    inputBinding:
      position: 998
      shellQuote: false
      prefix: '&& gzip *wig'
    default: true

  outputfolder:
    type: string?
    label: "Output directory name"
    inputBinding:
      position: 999
      shellQuote: false
      separate: false
      prefix: ' && nameoffolder="'
    default: "SICER_out"

  verifymove:
    type: boolean?
    label: "Move files to new directory"
    inputBinding:
      position: 1000 
      shellQuote: false
      prefix: '" && mkdir -p $nameoffolder && mv *W200* $nameoffolder' 
    default: true


outputs:
  sicerDir:
    type: Directory
    label: "Output directory"
    outputBinding:
      glob: |
        ${
          return inputs.outputfolder;
        }


