#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [igvtools, toTDF]
#[java, -jar, /research/rgs01/project_space/abrahgrp/Software_Dev_Sandbox/common/madetunj/software/igvtools/igvtools_2.3.2.jar, toTDF]
class: CommandLineTool
label: IGVTOOLS - convert to TDF
doc: |
  java -jar igvtools/igvtools_2.3.2.jar toTDF <wig> <tdf> hg19
  igvtools toTDF <wig> <tdf> hg19


hints:
  DockerRequirement:
    dockerPull: madetunj/igvtools:v2.8.2


requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var var_output_name = function() {
      if (inputs.wigfile != null) {
         return inputs.wigfile.basename.split('wig')[0]+'tdf';
      }
   };


inputs:
  wigfile:
    label: "Wig file"
    type: File
    inputBinding:
      position: 1
  
  totdffile:
    label: "TDF file name"
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
    label: "Genome name"
    type: string?
    default: "hg19"
    inputBinding:
      position: 3


outputs:
  outfile:
    label: "Output directory"
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.totdffile == "") {
            return var_output_name();
          } else {
            return inputs.totdffile;
          } 
        }
