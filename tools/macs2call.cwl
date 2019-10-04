#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [/hpcf/apps/python/install/2.7.5/bin/macs2, callpeak]
class: CommandLineTool

label: MACS - Model based Analysis from ChiP-Seq
doc: |
  macs2 callpeak -t $file.sorted.bam -g $specie -f BAM -p 1e-9 --keep-dup=auto -n $file\_p9_kd-auto

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

  expressionLib:
  - var var_output_name = function() {
      return inputs.treatmentfile.nameroot.split('.').slice(0,-1).join('.')+'_p9_kd-'+inputs.keep_dup+'_nm';
   };

inputs:
  treatmentfile:
    type: File
    inputBinding:
      prefix: -t
      position: 1
  
  species:
    type: string?
    default: hs
    inputBinding:
      prefix: -g

  fileformat:
    type: string?
    default: BAM
    inputBinding:
      prefix: -f
  
  pvalue:
    type: string?
    default: '1e-9'
    inputBinding:
      prefix: -p

  keep_dup:
    type: string?
    default: auto
    inputBinding:
      prefix: --keep-dup

  bedgraph:
    type: boolean?
    default: true
    inputBinding: 
      prefix: --bdg

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
      glob: '*narrowPeak'

  peaksxlsfile:
    type: File
    outputBinding:
      glob: '*peaks.xls'

  summitsfile:
    type: File
    outputBinding:
      glob: '*summits.bed'

  bedgraphfile:
    type: File
    outputBinding:
      glob: '*treat_pileup.bdg'


doc: |
  Model-based Analysis for ChIP-Sequencing
    Usage: macs2call.cwl [-h] [--pvalue PVALUE] [--keep-dup KEEP_DUPLICATES] [--outputname OUTPUTNAME]

    Options: --pvalue		STRING 	P-value
             --outputname	STRING	experiment outfile file name
             --keep-dup		STRING	Keep duplicates (default: auto)
             --treatmentfile	FILE	input CHIP bam file
                
 
$namespaces:
  s: http://schema.org/
 
 
$schemas:
 - https://schema.org/docs/schema_org_rdfa.html
