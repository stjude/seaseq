#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [bamToBed]
class: CommandLineTool

label: convert bam to bed
doc: |
  bamToBed -i <bam file> > <bed file>

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

  expressionLib:
  - var var_output_name = function() {
      if (inputs.infile != null) {
         return inputs.infile.nameroot.split('.bam')[0]+'.bam2bed.bed';
      }
   };

inputs:
  infile:
    type: File
    inputBinding:
      prefix: '-i'
      position: 1
  
  outfile:
    type: string?
    default: ""

stdout: |
  ${
    if (inputs.outfile == "") {
      return var_output_name();
    } else {
      return inputs.outfile;
    }
  }

outputs:
  outfile:
    type: stdout
