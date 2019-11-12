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
- class: InitialWorkDirRequirement
  listing: [ $(inputs.peaksbedfile) ]

inputs:
  inputdir:
    type: string?
    default: "./"
    inputBinding:
      position: 1
  
  peaksbedfile:
    type: File
    inputBinding:
      position: 2
      valueFrom: $(self.basename)

  outputdir:
    type: string?
    default: "./"
    inputBinding:
      position: 3

  species: 
    type: string?
    default: "hg19"
    inputBinding:
      position: 4
  
  redundancy:
    type: int?
    default: 1
    inputBinding:
      position: 5 

  window:
    type: int?
    default: 200
    inputBinding:
      position: 6

  fragment_size:
    type: int?
    default: 150
    inputBinding:
      position: 7

  genome_fraction:
    type: double?
    default: 0.86
    inputBinding:
      position: 8

  gapsize:
    type: int?
    default: 200
    inputBinding:
      position: 9

  evalue:
    type: int?
    default: 100
    inputBinding:
      position: 10

outputs:
  islandbed:
    type: File
    outputBinding:
      glob: '*islandfiltered.bed'

  scoreisland:
    type: File
    outputBinding:
      glob: '*scoreisland'

  graph:
    type: File
    outputBinding:
      glob: '*graph'
