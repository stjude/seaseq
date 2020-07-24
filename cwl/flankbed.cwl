#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [ flanking.pl ]
class: CommandLineTool
doc: |
  flanking.pl summits.bed 50 > summits-flank50.bed


hints:
  DockerRequirement:
    dockerPull: madetunj/seaseq:v0.0.1


requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var var_output_name = function() {
      if (inputs.bedfile != null) {
         return inputs.bedfile.nameroot+'-flank'+inputs.flank+'.bed';
      }
   };


inputs:
  flank:
    type: int?
    label: "bp distance from summit"
    default: 50
    inputBinding:
      prefix: -f

  bedfile:
    type: File
    label: "Summit BED file"
    inputBinding:
      prefix: -i
  
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
    label: "Output BED file"
