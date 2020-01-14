#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [intersectBed, -v]
class: CommandLineTool

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
    inputBinding:
      prefix: '-a'
      position: 1
  
  blacklistfile:
    type: File
    inputBinding:
      prefix: '-b'
      position: 2


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
