#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [java, -jar, igvtools/igvtools_2.3.2.jar, toTDF]
class: CommandLineTool

label: IGVTOOLS - convert to TDF
doc: |
  java -jar igvtools/igvtools_2.3.2.jar toTDF <wig[not gziped]> <tdf> hg19

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

  expressionLib:
  - var var_output_name = function() {
      if (inputs.infile != null) {
         return inputs.wigfile.basename.split('wig')[0]+'tdf';
      }
   };

inputs:
  wigfile:
    type: File
    inputBinding:
      position: 1
  
  totdffile:
    type: string
    inputBinding:
      position: 2
      valueFrom: |
        ${
            if (self == ""){
              return var_output_name();
            } else {
              return self;
            }
        }
    default: ""

  genome:
    type: string?
    default: "hg19"
    inputBinding:
      position: 3

outputs:
  outfile:
    type: File
    outputBinding:
      glob: '*.tdf'
