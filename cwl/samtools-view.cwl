#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [samtools, view]
class: CommandLineTool
label: convert sam to bam file


hints:
  DockerRequirement:
    dockerPull: madetunj/samtools:v1.9 


requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var var_output_name = function() {
      if (inputs.infile != null) {
         return inputs.infile.nameroot+'.bam.bam';
      }
   };


inputs:
  infile:
    type: File
    label: "SAM file"
    inputBinding:
      prefix: -b

  outputfile:
    label: "Output file name"
    type: string?
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
    label: "BAM file"

