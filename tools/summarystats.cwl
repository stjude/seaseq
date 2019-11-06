#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [ summaryfacts.pl ]
class: CommandLineTool
label: QC on the peaks file
doc: perl summaryfacts.pl -b <bam file> -p <peaks bed> -bamflag <bamflagstat> -rmdupflag <rmdupflagstat> -bkflag <bklistflagstat> -fqc <fastqczipfile> -outfile <outputfilename>

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var var_output_name = function() {
      if (inputs.outfile == ""){
        if (inputs.peaksbed != null) { return inputs.peaksbed.nameroot+'.txt'; }
        else if (inputs.bamflag != null) { return inputs.bamflag.nameroot+'.txt'; }
        else if (inputs.rmdupflag != null) { return inputs.rmdupflag.nameroot+'.txt'; }
        else if (inputs.bkflag != null) { return inputs.bkflag.nameroot+'.txt'; }
      }
   };

- class: InitialWorkDirRequirement
  listing: [ $(inputs.bamfile) ]

inputs:
  bamfile:
    type: File?
    inputBinding:
      prefix: -b
      valueFrom: $(self.basename)

  peaksbed:
    type: File?
    inputBinding:
      prefix: -p

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
