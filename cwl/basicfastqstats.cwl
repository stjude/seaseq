#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [basicfastqstats.sh]
class: CommandLineTool


hints:
  DockerRequirement:
    dockerPull: madetunj/seaseq:v0.0.1


requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var var_output_name = function() {
      if (inputs.fastqfile != null) {
         return inputs.fastqfile.nameroot.split('.fastq')[0]+'-fastq.metrics.txt';
      }
   };


inputs:
  fastqfile:
    type: File
    label: "FastQ file"
    inputBinding:
      position: 1

  outfile:
    type: string?
    label: "Output file name"
    inputBinding:
      position: 2
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
  metrics_out:
    type: File
    label: "Output file"
    outputBinding:
      glob: |
        ${
          if (inputs.outfile == "") {
            return var_output_name();
          } else {
            return inputs.outfile;
          }
        }
