#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [ summaryfactsv2.pl ]
class: CommandLineTool
label: QC on the peaks file
doc: perl summaryfacts.pl -b <bam file> -p <peaks bed> -bamflag <bamflagstat> -rmdupflag <rmdupflagstat> -bkflag <bklistflagstat> -fqc <fastqczipfile> -outfile <outputfilename>

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var var_output_name = function() {
      if (inputs.outfile == ""){
        if (inputs.bamflag != null) { return inputs.bamflag.nameroot+'-stats.out'; }
        else if (inputs.rmdupflag != null) { return inputs.rmdupflag.nameroot+'-stats.out'; }
        else if (inputs.bkflag != null) { return inputs.bkflag.nameroot+'-stats.out'; }
      }
   };

inputs:
  bambed:
    type: File?
    inputBinding:
      prefix: -b

  sppfile:
    type: File?
    inputBinding:
      prefix: -s

  countsfile:
    type: File?
    inputBinding:
      prefix: -c

  peaksxls:
    type: File?
    inputBinding:
      prefix: -px

  bamflag:
    type: File?
    inputBinding:
      prefix: '-bamflag'

  rmdupflag:
    type: File?
    inputBinding:
      prefix: '-rmdupflag'

  bkflag:
    type: File?
    inputBinding:
      prefix: '-bkflag'

  fastqczip:
    type: File?
    inputBinding:
      prefix: '-fqc'

  fastqmetrics:
    type: File?
    inputBinding:
      prefix: '-fx'

  rosedir:
    type: Directory?
    inputBinding:
      prefix: '-rose'

  outfile:
    type: string?
    inputBinding:
      prefix: '-outfile'
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
  statsfile:
    type: File
    outputBinding: 
      glob: '*stats.out'

  htmlfile:
    type: File
    outputBinding:
      glob: '*stats.html'

  textfile:
    type: File
    outputBinding:
      glob: '*stats.txt'
