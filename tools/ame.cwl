#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: ame
class: CommandLineTool

label: AME - Analysis of Motif Enrichment
doc: |
  ame <convert fasta> <motif-databases>

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  
  expressionLib:
  - var var_output_name = function() {
      return inputs.convertfasta.nameroot+'-ame';
    };

inputs:
  convertfasta:
    type: File
    inputBinding:
      position: 1
  
  motifdatabases:
    type: 
      type: array
      items: string
    inputBinding:
      position: 2

  outfile_txt:
    type: string?
    inputBinding:
      position: 1000
      shellQuote: false
      prefix: '&& mv ame_out/ame.txt'
      valueFrom: |
        ${
            if (self == ""){
              return var_output_name()+'.txt';
            } else {
              return self;
            }
        }
    default: ""

  outfile_html:
    type: string?
    inputBinding:
      position: 999
      shellQuote: false
      prefix: '&& mv ame_out/ame.html'
      valueFrom: |
        ${
            if (self == ""){
              return var_output_name()+'.html';
            } else {
              return self;
            }
        }
    default: ""


outputs:
  outfile:
    type: File
    outputBinding:
      glob: '*ame.txt'

  htmlfile:
    type: File
    outputBinding:
      glob: '*ame.html'
