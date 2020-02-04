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
  - var var_output_prefix = function() {
      if (inputs.outputfolder == "") {
        if (inputs.outputname == "") {
          return inputs.treatmentfile.basename.split('.').slice(0,-1).join('.')+'-nm'+'/'+inputs.treatmentfile.basename.split('.').slice(0,-1).join('.')+'_nm';
        }
        else {
          return inputs.treatmentfile.basename.split('.').slice(0,-1).join('.')+'-nm'+'/'+inputs.outputname;
        }
      }
      else {
        if (inputs.outputname == "") {
          return inputs.outputfolder+'/'+inputs.treatmentfile.basename.split('.').slice(0,-1).join('.')+'_nm';
        }
        else {
          return inputs.outputfolder+'/'+inputs.outputname;
        }
      }
   };
  
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

  outputfolder:
    type: string?
    inputBinding:
      position: 999
      shellQuote: false
      separate: false
      prefix: ' && nameoffolder="'
      valueFrom: |
        ${
            if (self == ""){
              return var_output_folder();
            } else {
              return self;
            }
        }
    default: ""

  verifymove:
    type: boolean?
    inputBinding:
      position: 1000 
      shellQuote: false
      prefix: '" && mkdir -p $nameoffolder && mv *_nm* $nameoffolder' 
    default: true


outputs:
  peaksbedfile:
    type: File
    outputBinding:
      glob: |
        ${
          return var_output_prefix() + '_peaks.bed';
        }

  peaksxlsfile:
    type: File
    outputBinding:
      glob: |
        ${
          return var_output_prefix() + '_peaks.xls';
        }

  summitsfile:
    type: File
    outputBinding:
      glob: |
        ${
          return var_output_prefix() + '_summits.bed';
        }

  wigfile:
    type: File
    outputBinding:
      glob: |
        ${
          return var_output_prefix() + '_MACS_wiggle/treat/*_treat_afterfiting_all.wig.gz';
        }

  macsDir:
    type: Directory
    outputBinding:
      glob: |
        ${
          if (inputs.outputfolder == "") {
            return var_output_folder();
          } else {
            return inputs.outputfolder;
          }
        }
