#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [samtools, index]
class: CommandLineTool
label: index bam file


hints:
  DockerRequirement:
    dockerPull: madetunj/samtools:v1.9 


requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.infile) ]


inputs:
  infile:
    type: File
    label: "BAM file"
    inputBinding:
      position: 1
      valueFrom: $(self.basename)


outputs:
  outfile:
    type: File
    label: "BAM index output file name"
    secondaryFiles: .bai
    outputBinding:
     glob: $(inputs.infile.basename)

