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
         return inputs.infile.nameroot.split('.rmdup')[0]+'.rm-blist.bam';
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
    type: string
    inputBinding:
      position: 3
      prefix: '>'
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
      glob: '*rm-blist.bam'
