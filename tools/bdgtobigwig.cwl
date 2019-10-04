#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [ bedGraphToBigWig ]
class: CommandLineTool
#module load ucsc/041619

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

  expressionLib:
  - var var_output_name = function() {
      if (inputs.infile != null) {
         return inputs.infile.basename.split('_nm')[0]+'.bw';
      }
   };

inputs:
  infile:
    label: peak calling file in bedgraph format
    type: File
    inputBinding:
      position: 1

  chromsizes:
    label: chromosome sizes in a two-column file
    type: File
    inputBinding:
      position: 2

  outfile:
    label: outputfilename
    type: string?
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
    outputBinding: 
      glob: '*bw'

