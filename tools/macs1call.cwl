#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [macs14]
class: CommandLineTool

label: MACS - Model based Analysis from ChiP-Seq
doc: |
  macs14 -t $file.bam -w -S --space=50 -p 1e-9 --keep-dup=auto -n $file\_p9_kd-auto

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var var_output_name = function() {
      return inputs.treatmentfile.basename.split('.').slice(0,-1).join('.')+'_p9_kd-'+inputs.keep_dup;
   };
  - var var_output_folder = function() {
      return inputs.treatmentfile.basename.split('.').slice(0,-1).join('.')+'-p9_kd-'+inputs.keep_dup;
   };
- class: InitialWorkDirRequirement
  listing: [ $(inputs.treatmentfile) ]


inputs:
  treatmentfile:
    type: File
    inputBinding:
      prefix: -t
      position: 1

  pvalue:
    type: string?
    default: '1e-9'
    inputBinding:
      prefix: -p

  keep_dup:
    type: string?
    default: auto
    inputBinding:
      prefix: --keep-dup=
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

  space:
    type: int?
    default: 50
    inputBinding:
      prefix: --space=
      separate: false

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
      prefix: '&& mv *_p9_kd-* ' 
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
    Usage: macs1call.cwl [-h] --treatmentfile TREATMENTFILE [--pvalue PVALUE] 
                    [--keep-dup KEEP_DUPLICATES] [--outputname OUTPUTNAME]

    Options: --pvalue		STRING 	P-value
             --outputname	STRING	experiment outfile file name
             --keep-dup		STRING	Keep duplicates (default: auto)
             --treatmentfile	FILE	input CHIP bam file
                
 
$namespaces:
  s: http://schema.org/
 
 
$schemas:
 - https://schema.org/docs/schema_org_rdfa.html
