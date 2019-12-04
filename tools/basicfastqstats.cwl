#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [basicfastqstats.sh]
class: CommandLineTool

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
    label: "FastQfiles"
    inputBinding:
      position: 1

  outputresults:
    type: string?
    label: "Output File"
    inputBinding:
      position: 2
      prefix: '>'
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
    outputBinding:
      glob: '*metrics.txt'
