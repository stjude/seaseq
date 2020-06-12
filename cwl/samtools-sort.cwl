#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [samtools, sort]
class: CommandLineTool
label: sort bam file

hints:
  DockerRequirement:
    dockerPull: madetunj/samtools:v1.9

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

  expressionLib:
  - var var_output_name = function() {
      if (inputs.infile != null) {
         return inputs.infile.nameroot.split('.bam')[0]+'.sorted.bam';
      }
   };

inputs:
  infile:
    label: "BAM file"
    type: File
    inputBinding:
      position: 1

  outputfile:
    type: string
    label: "Output file name"
    inputBinding:
      prefix: '-o'
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
    label: "Sorted BAM output file"
    outputBinding:
      glob: |
        ${
          if (inputs.outputfile == "") {
            return var_output_name();
          } else {
            return inputs.outputfile;
          } 
        }
