#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [ bedtools, getfasta ]
class: CommandLineTool

doc: |
  bedtools getfasta -fi <peads bed> -fo <fasta outputfile>


requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

  expressionLib:
  - var var_output_name = function() {
      if (inputs.bedfile != null) {
         return inputs.bedfile.nameroot+'.fa';
      }
   };

inputs:
  reference:
    type: Directory
    inputBinding:
      prefix: -fi
      position: 1
      valueFrom: |
        ${
            for (var i = 0; i < self.listing.length; i++) {
                if (self.listing[i].path.split('.').slice(-1) == 'fa') {
                  return self.listing[i].path;
                }
            }
            return null;
        }

  bedfile:
    type: File
    inputBinding:
      position: 2
      prefix: -bed
  
  outfile:
    type: string
    inputBinding:
      position: 3
      prefix: -fo
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
      glob: '*fa'
