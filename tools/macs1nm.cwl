#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [macs14]
class: CommandLineTool

label: MACS - Model based Analysis from ChiP-Seq
doc: |
  macs14 -t $file.bam -w -S --shiftsize=100 --space=50 --nomodel -n $file\_nm

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var var_output_name = function() {
      return inputs.treatmentfile.basename.split('.').slice(0,-1).join('.')+'_nm';
   };
  - var var_output_folder = function() {
      return inputs.treatmentfile.basename.split('.').slice(0,-1).join('.')+'-nm';
   };
- class: InitialWorkDirRequirement 
  listing: [ $(inputs.treatmentfile) ]


inputs:
  treatmentfile:
    type: File
    inputBinding:
      prefix: -t
      position: 1

  space:
    type: int?
    default: 50
    inputBinding:
      prefix: --space=
      separate: false

  nomodel:
    type: boolean?
    default: true
    inputBinding:
      prefix: --nomodel

  shiftsize:
    type: int?
    default: 200
    inputBinding:
      prefix: --shiftsize=
      separate: false
  
  wiggle:
    type: boolean?
    default: true
    inputBinding:
      prefix: -w
  
  single_profile:
    type: boolean?
    default: true
    inputBinding:
      prefix: -S

  outputname:
    type: string?
    inputBinding:
      prefix: -n
      valueFrom: |
        ${
            if (self == ""){
              return var_output_name();
            } else {
              return self;
            }
        }
    default: ""

  movefiles:
    type: string?
    inputBinding:
      position: 999
      shellQuote: false
      prefix: '&& mkdir -p '
      valueFrom: |
        ${
            if (self == ""){
              return var_output_folder();
            } else {
              return self;
            }
        }
    default: ""

  movefiles2:
    type: string?
    inputBinding:
      position: 1000
      shellQuote: false
      prefix: '&& mv *_nm* '
      valueFrom: |
        ${
            if (self == ""){
              return var_output_folder();
            } else {
              return self;
            }
        }
    default: ""


outputs:
  peaksbedfile:
    type: File
    outputBinding:
      glob: '*peaks.bed'

  peaksxlsfile:
    type: File
    outputBinding:
      glob: '*peaks.xls'

  summitsfile:
    type: File
    outputBinding:
      glob: '*summits.bed'

  wigfile:
    type: File
    outputBinding:
      glob: '*_treat_afterfiting_all.wig.gz'

  macsDir:
    type: Directory
    outputBinding:
      glob: |
        ${
          if (inputs.movefiles == "") {
            return var_output_folder();
          } else {
            return inputs.movefiles;
          }
        }


doc: |
  Model-based Analysis for ChIP-Sequencing
    Usage: macs1nm.cwl [-h] [--outputname OUTPUTNAME]

    Options: --outputname	STRING	experiment outfile file name
             --treatmentfile	FILE	input CHIP bam file
                
 
$namespaces:
  s: http://schema.org/
 
 
$schemas:
 - https://schema.org/docs/schema_org_rdfa.html
