#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: meme-chip
class: CommandLineTool
label: MEME-ChIP performs comprehensive motif analysis (including motif discovery) 
doc: |
  meme-chip <convert fasta>

  
hints:
  DockerRequirement:
    dockerPull: madetunj/memesuite:v5.1.1


requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var var_output_name = function() {
      return 'bklist'+inputs.convertfasta.nameroot.split('bklist').slice(-1)+'-memechip_out';
   };


inputs:
  convertfasta:
    label: "BED converted FASTA file"
    type: File
    inputBinding:
      position: 1000

  spamo-skip:
    label: "remove SPAced MOtif analysis"
    type: boolean?
    default: false
    inputBinding:
      position: 1
      prefix: '-spamo-skip'

  fimo-skip:
    label: "remove Find Individual Motif Occurences"
    type: boolean?
    default: false
    inputBinding:
      position: 2
      prefix: '-fimo-skip'
    
  outputdir:
    label: "MEME-chip output directory name"
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
