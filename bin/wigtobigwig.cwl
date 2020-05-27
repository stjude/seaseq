#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [wigToBigWig, -clip]
class: CommandLineTool
label: convert ascii format wig file to binary big wig format

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

  expressionLib:
  - var var_output_name = function() {
      if (inputs.infile != null) {
         return inputs.infile.nameroot.split('wig')[0]+'bw';
      }
   };

inputs:
  infile:
    type: File
    label: "Wig file"
    inputBinding:
      position: 1

  chromsizes:
    type: File
    label: "Chromosome sizes tab file"
    inputBinding:
      position: 2

  outputfile:
    type: string
    label: "BigWig output file name"
    inputBinding:
      position: 3
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
    label: "Output file"
    outputBinding:
      glob: |
        ${
          if (inputs.outputfile == "") {
            return var_output_name();
          } else {
            return inputs.outputfile;
          } 
        }
