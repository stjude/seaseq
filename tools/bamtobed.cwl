#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [bamToBed]
class: CommandLineTool

label: convert bam to bed
doc: |
  bamToBed -i <bam file> > <bed file>

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

  expressionLib:
  - var var_output_name = function() {
      if (inputs.infile != null) {
         return inputs.infile.nameroot.split('.bam')[0]+'.bam2bed.bed';
      }
   };

inputs:
  infile:
    label: "BAM file"
    type: File
    inputBinding:
      prefix: '-i'
      position: 1
  
  outputfile:
    type: string?
    label: "Output BED file name"
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
    label: "Output BED file"
