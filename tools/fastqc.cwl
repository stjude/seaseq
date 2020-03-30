
#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [fastqc]
class: CommandLineTool
label: QC on reads or bam file

requirements:
- class: InitialWorkDirRequirement
  listing: [ $(inputs.infile) ]

inputs:
  infile:
    type: File
    inputBinding:
      position: 1
  
  outputdirectory:
    type: string?
    inputBinding:
      position: 2
      prefix: '-o'
    default: './'

outputs:
  htmlfile:
    type: File
    outputBinding: 
      glob: '*_fastqc.html'

  zipfile:
    type: File
    outputBinding:
      glob: '*_fastqc.zip'
