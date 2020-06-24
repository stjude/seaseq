#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [samtools, flagstat]
class: CommandLineTool
label: SamTools flagstat
doc: |
  samtools flagstat $BAM > flagstat.txt


hints:
  DockerRequirement:
    dockerPull: madetunj/samtools:v1.9


requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var var_output_name = function() {
      return inputs.infile.nameroot+'-flagstat.txt';
   };


inputs:
  infile:
    type: File
    label: "BAM file"
    inputBinding:
      position: 1
  
  outputfile:
    type: string?
    label: "Output file name"
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
    label: "FlagStats standard output"
