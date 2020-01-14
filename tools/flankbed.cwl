#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [ flanking.pl ]
class: CommandLineTool

doc: |
  flanking.pl summits.bed 50 > summits-flank50.bed


requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

  expressionLib:
  - var var_output_name = function() {
      if (inputs.bedfile != null) {
         return inputs.bedfile.nameroot+'-flank'+inputs.flank+'.bed';
      }
   };

inputs:
  flank:
    type: int?
    default: 50
    inputBinding:
      prefix: -f

  bedfile:
    type: File
    inputBinding:
      prefix: -i
  
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
