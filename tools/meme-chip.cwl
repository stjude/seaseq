#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: meme-chip
class: CommandLineTool

label: MEME-ChIP performs comprehensive motif analysis (including motif discovery) 
doc: |
  meme-chip <convert fasta> 

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

  expressionLib:
  - var var_output_name = function() {
      return inputs.convertfasta.nameroot+'-meme';
    };

inputs:
  convertfasta:
    type: File
    inputBinding:
      position: 1

  outfile_txt:
    type: string?
    inputBinding:
      position: 1000
      shellQuote: false
      prefix: '&& mv memechip_out/meme_out/meme.txt'
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
      prefix: '&& mv memechip_out/meme_out/meme.html'
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
      glob: '*meme.txt'

  htmlfile:
    type: File
    outputBinding:
      glob: '*meme.html'
