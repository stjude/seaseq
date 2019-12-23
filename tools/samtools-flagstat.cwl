#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [samtools, flagstat]
class: CommandLineTool

label: SamTools flagstat
doc: |
  samtools flagstat $BAM > flagstat.txt

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
    inputBinding:
      position: 1
  
  outfile:
    type: string
    inputBinding:
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
  outfile:
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
