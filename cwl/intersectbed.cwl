#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [intersectBed]
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: madetunj/bedtools:v2.25.0

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
    label: "Peaks BED file"
    inputBinding:
      prefix: '-a'
      position: 1

  bamtobed:
    type: File
    label: "Peaks BAM to BED file"
    inputBinding:
      prefix: '-b'
      position: 1
  
  countoverlap:
    type: boolean?
    label: "Number of overlap with bamtobed file"
    inputBinding:
      prefix: '-c'
      position: 1
    default: true

  sorted:
    type: boolean?
    label: "Sort bed files"
    inputBinding: 
      prefix: '-sorted'
      position: 1
    default: true

  outputfile:
    type: string?
    label: "Output file name"
    default: ""

stdout: |
  ${
    if (inputs.outputfile == "") {
      return var_output_name();
    } else {
      return inputs.outputfile;
    }
  }

outputs:
  outfile:
    type: stdout
    label: "Output file"
