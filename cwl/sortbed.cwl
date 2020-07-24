#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [sort-bed]
class: CommandLineTool
label: Using bedops to sort bed file
doc: |
  sort-bed <bed file> > <sorted bed file>


hints:
  DockerRequirement:
    dockerPull: madetunj/bedops:v2.4.37


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
    label: "BED file"
    inputBinding:
      position: 1
  
  outputfile:
    type: string?
    label: "sorted BED output file name"
    default: ""


stdout: |
  ${
    if (inputs.outputfile == "") {
      return var_output_name();
    } else {
      return inputs.outputfile;
    }
  }


outputs:
  outfile:
    type: stdout
    label: "sorted BED file"
