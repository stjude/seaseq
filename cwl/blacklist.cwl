#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [intersectBed, -v]
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: madetunj/bedtools:v2.25.0

label: bedtools intersect, to remove blacklist
doc: |
  intersectBed -v -a KOPTK1_DMSO.rmdup.bam -b ~/.genomes/hg19/hg19-blacklist.v2.bed > ooo

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

  expressionLib:
  - var var_output_name = function() {
      if (inputs.infile != null) {
         return inputs.infile.nameroot.split('.bam')[0]+'.bklist.bam';
      }
   };

inputs:
  infile:
    type: File
    label: "BAM file"
    inputBinding:
      prefix: '-a'
      position: 1
  
  blacklistfile:
    type: File
    label: "Blacklist file"
    inputBinding:
      prefix: '-b'
      position: 2


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
    label: "Output BAM file"
