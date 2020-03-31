#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [run_spp.R]
class: CommandLineTool
label: Quality metrics using PhantomPeaksQual tool
doc: |
  run_spp.R -c=<bam> -savp -out=<outfile name> 

requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var var_output_name = function() {
      return inputs.infile.nameroot.split('.bam')[0]+'-spp.out';
   };

inputs:
  infile:
    type: File
    inputBinding:
      prefix: '-c='
      separate: false
      position: 1

  crosscorr:
    type: boolean?
    inputBinding:
      position: 2
      prefix: '-savp'
    default: true

  outfile:
    type: string?
    default: ""
    inputBinding:
      prefix: '-out='
      separate: false
      position: 3
      valueFrom: |
        ${
            if (self == ""){
              return var_output_name();
            } else {
              return self;
            }
        }

outputs:
  spp_out:
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
