#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [intersectBed]
class: CommandLineTool

label: number of overlap of A with B
doc: |
  intersectBed -sorted -a <peaksbed> -b <sorted bamtobed> -c > <output countfile>

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

  expressionLib:
  - var var_output_name = function() {
      if (inputs.peaksbed != null) {
         return inputs.peaksbed.nameroot.split('.bed')[0]+'-counts.txt';
      }
   };

inputs:
  peaksbed:
    type: File
    inputBinding:
      prefix: '-a'
      position: 1

  bamtobed:
    type: File
    inputBinding:
      prefix: '-b'
      position: 1
  
  countoverlap:
    type: boolean?
    inputBinding:
      prefix: '-c'
      position: 1
    default: true

  sorted:
    type: boolean?
    inputBinding: 
      prefix: '-sorted'
      position: 1
    default: true

  outfile:
    type: string?
    default: ""

stdout: |
  ${
    if (inputs.outfile == "") {
      return var_output_name();
    } else {
      return inputs.outfile;
    }
  }

outputs:
  outfile:
    type: stdout
