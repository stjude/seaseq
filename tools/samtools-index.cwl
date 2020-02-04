#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [samtools, index]
class: CommandLineTool
label: index bam file

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.infile) ]

inputs:
  infile:
    type: File
    inputBinding:
      position: 1
      valueFrom: $(self.basename)

outputs:
  outfile:
    type: File
    secondaryFiles: .bai
    outputBinding:
     glob: $(inputs.infile.basename)

