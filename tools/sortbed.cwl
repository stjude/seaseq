#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [sort-bed]
class: CommandLineTool

label: Using bedops to sort bed file
doc: |
  sort-bed <bed file> > <sorted bed file>

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

  expressionLib:
  - var var_output_name = function() {
      if (inputs.infile != null) {
         return inputs.infile.nameroot.split('.bed')[0]+'.sorted.bed';
      }
   };

inputs:
  infile:
    type: File
    inputBinding:
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
