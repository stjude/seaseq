#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [wigToBigWig, -clip]
class: CommandLineTool
label: convert ascii format wig file to binary big wig format

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

  expressionLib:
  - var var_output_name = function() {
      if (inputs.infile != null) {
         return inputs.infile.nameroot.split('wig')[0]+'bw';
      }
   };

inputs:
  infile:
    type: File
    inputBinding:
      position: 1

  chromsizes:
    type: File
    inputBinding:
      position: 2

  outfile:
    type: string
    inputBinding:
      position: 3
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
      glob: |
        ${
          if (inputs.outfile == "") {
            return var_output_name();
          } else {
            return inputs.outfile;
          } 
        }


doc: |
  Convert ascii format wig file to binary big wig format
    Usage: wigtobigwig.cwl [-h] --infile INFILE --chromsizes CHROMSIZES --outfile OUTFILE

    Options: --infile		FILE 	peak calling file in wiggle format
             --chromsizes	FILE	chromosome sizes in a two-column file
             --outfile		INT	outputfilename (optional)
                
 
$namespaces:
  s: http://schema.org/
 
 
$schemas:
 - https://schema.org/docs/schema_org_rdfa.html
