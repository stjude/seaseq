#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: macs14
class: CommandLineTool

label: MACS - Model based Analysis from ChiP-Seq
doc: |
  bsub -K -R \"select[rhel7]\" macs14 -t $file.sorted.bam -w -S --space=50 --nomodel --shiftsize=200 -n $file\_nm" >> $outfile

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

  expressionLib:
  - var var_output_name = function() {
      return inputs.treatmentfile.basename.split('.').slice(0,-1).join('.')+'_nm';
   };

inputs:
  treatmentfile:
    type: File
    inputBinding:
      prefix: -t
      position: 1
  
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

  outputname:
    type: string
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
    Usage: macs.cwl [-h] [--nomodel] --treatmentfile TREATMENTFILE --outputname OUTPUTNAME 
                [--shiftsize SHIFTSIZE] [--single_profile] [--space SPACE] [--nomodel] 
                [--wiggle]

    Options: --treatmentfile	FILE 	input CHiP bam file
             --outputname	FILE	experiment outfile file name 
             --shiftsize	INT	bp shift size (default: 200)
             --space		INT	saved wiggle file resolution (default: 50)
             --wiggle		BOOLEAN	save extended fragment pileup into a wiggle file (default: true)
             --nomodel		BOOLEAN	build shifting model (default: true)
             --single_profile	BOOLEAN	entire genome for treatment and input (default: true)
                
 
$namespaces:
  s: http://schema.org/
 
 
$schemas:
 - https://schema.org/docs/schema_org_rdfa.html
