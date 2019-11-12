#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [macs14]
class: CommandLineTool

label: MACS - Model based Analysis from ChiP-Seq
doc: |
  macs14 -t $file.bam -w -p 1e-9 --keep-dup=auto -n $file\_p9_kd-auto

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

  expressionLib:
  - var var_output_name = function() {
      return inputs.treatmentfile.basename.split('.').slice(0,-1).join('.')+'_p9_kd-'+inputs.keep_dup+'_nm';
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
