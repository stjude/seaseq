#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: ame
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: madetunj/memesuite:v5.1.1 

label: AME - Analysis of Motif Enrichment
doc: |
  ame <convert fasta> <motif-databases>

requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var var_output_name = function() {
      return 'bklist'+inputs.convertfasta.nameroot.split('bklist').slice(-1)+'-ame_out';
   };

inputs:
  convertfasta:
    label: "BED converted FASTA file"
    type: File
    inputBinding:
      position: 999
  
  motifdatabases:
    label: "MEME motif databases to identify motif enrichment"
    type: 
      type: array
      items: File
    inputBinding:
      position: 1000

  outputdir:
    label: "AME output directory name"
    type: string?
    inputBinding:
      position: 1
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
    label: "Output directory"
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

