#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [samtools, sort]
class: CommandLineTool
label: sort bam file

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

  expressionLib:
  - var var_output_name = function() {
      if (inputs.infile != null) {
         return inputs.infile.nameroot.split('.bam')[0]+'.sorted.bam';
      }
   };

inputs:
  infile:
    type: File
    inputBinding:
      position: 1

  outfile:
    type: string
    inputBinding:
      prefix: '-o'
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
      glob: |
        ${
          if (inputs.outfile == "") {
            return var_output_name();
          } else {
            return inputs.outfile;
          } 
        }
