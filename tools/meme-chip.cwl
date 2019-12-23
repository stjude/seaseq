#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: meme-chip
class: CommandLineTool

label: MEME-ChIP performs comprehensive motif analysis (including motif discovery) 
doc: |
  meme-chip <convert fasta>
  
requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var var_output_name = function() {
      return inputs.convertfasta.nameroot+'-memechip_out';
   };

inputs:
  convertfasta:
    type: File
    inputBinding:
      position: 1000

  spamo-skip:
    type: boolean?
    default: false
    inputBinding:
      position: 1
      prefix: '-spamo-skip'

  fimo-skip:
    type: boolean?
    default: false
    inputBinding:
      position: 2
      prefix: '-fimo-skip'
    
  outputdir:
    type: string?
    inputBinding:
      position: 3
      prefix: -oc
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
  outDir:
    type: Directory
    outputBinding:
      glob: |
        ${
          if (inputs.outputdir == "") {
            return var_output_name();
          } else {
            return inputs.outputdir;
          } 
        }
