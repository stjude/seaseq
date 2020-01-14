#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [samtools, view]
class: CommandLineTool
label: convert sam to bam file

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

  expressionLib:
  - var var_output_name = function() {
      if (inputs.infile != null) {
         return inputs.infile.nameroot+'.bam.bam';
      }
   };

inputs:
  infile:
    type: File
    inputBinding:
      prefix: -b

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

