#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [run_spp.R]
class: CommandLineTool
label: Quality metrics using PhantomPeaksQual tool
doc: |
  run_spp.R -c=<bam> -savp -out=<outfile name> 

hints:
  DockerRequirement:
    dockerPull: madetunj/spp:v1.16.0

requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var var_output_name = function() {
      return inputs.infile.nameroot.split('.bam')[0]+'-spp.out';
   };
- class: InitialWorkDirRequirement
  listing: [ $(inputs.infile) ]

inputs:
  infile:
    type: File
    label: "BAM file"
    inputBinding:
      prefix: '-c='
      separate: false
      position: 1

  crosscorr:
    type: boolean?
    label: "save cross-correlation plot"
    inputBinding:
      position: 2
      prefix: '-savp'
    default: true

  outfile:
    type: string?
    label: "output file name"
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
    label: "output file"
    outputBinding:
      glob: |
        ${
          if (inputs.outfile == "") {
            return var_output_name();
          } else {
            return inputs.outfile;
          }
        }
