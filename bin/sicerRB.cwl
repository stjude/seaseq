#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [SICER-rb_custom.sh]
class: CommandLineTool

label: SICER - Broad Peaks
doc: |
  bsub -K -R \"select[rhel7]\" ~/.software/SICER_V1.1/SICER/SICER-rb.sh ./ <mapped converted bed> sicer hg19 1 200 150 0.86 0 100 >> $outfile
  [InputDir] [bed file] [OutputDir] [species] [redundancy threshold] [window size (bp)] [fragment size] [effective genome fraction] [gap size (bp)] [E-value]
  mkdir sicer5; ~/.software/SICER_V1.1/SICER/SICER-rb.sh ./ KOPTK1_DMSO.bam2bed.bed sicer5 hg19 1 200 150 0.86 200 100

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing: [ $(inputs.treatmentbedfile) ]

inputs:
  inputdir:
    type: string?
    label: "input files directory"
    default: "./"
    inputBinding:
      position: 1
  
  treatmentbedfile:
    type: File
    label: "BED file"
    inputBinding:
      position: 2
      valueFrom: $(self.basename)

  outputdir:
    type: string?
    label: "Output directory"
    default: "./"
    inputBinding:
      position: 3

  species: 
    type: string?
    label: "Genome name"
    default: "hg19"
    inputBinding:
      position: 4
  
  redundancy:
    type: int?
    label: "Number of identical reads allowed"
    default: 1
    inputBinding:
      position: 5 

  window:
    type: int?
    label: "Resolution of SICER"
    default: 200
    inputBinding:
      position: 6

  fragment_size:
    type: int?
    label: "Fragment shift size"
    default: 150
    inputBinding:
      position: 7

  genome_fraction:
    type: double?
    label: "Effective Genome Fraction"
    default: 0.86
    inputBinding:
      position: 8

  gapsize:
    type: int?
    label: "Gap size"
    default: 200
    inputBinding:
      position: 9

  evalue:
    type: int?
    label: "E-value"
    default: 100
    inputBinding:
      position: 10

  gzip_graph:
    type: boolean?
    label: "GZIP wig files"
    inputBinding:
      position: 998
      shellQuote: false
      prefix: '&& gzip *graph'
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
