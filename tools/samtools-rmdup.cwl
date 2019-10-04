#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [samtools, rmdup]
class: CommandLineTool
label: remove duplicates from bam file

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

  expressionLib:
  - var var_output_name = function() {
      if (inputs.infile != null) {
         return inputs.infile.nameroot.split('.sorted')[0]+'.rmdup.bam';
      }
   };

inputs:
  infile:
    type: File
    inputBinding:
      prefix: '-s'
      position: 1

  outfile:
    type: string
    inputBinding:
      position: 2
      valueFrom: |
        ${
            if (self == ""){
              return var_output_name();
            } else {
              return self;
            }
        }
    default: ""


outputs:
  outfile:
    type: File
    outputBinding: 
      glob: '*rmdup.bam'
