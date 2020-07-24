#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [macs14]
class: CommandLineTool
label: MACS1 - Model based Analysis from ChiP-Seq
doc: |
  macs14 -t $file.bam -w -S --space=50 -p 1e-9 --keep-dup=auto -n $file\_p9_kd-auto


hints:
  DockerRequirement:
    dockerPull: madetunj/macs:v1.4.2


requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var var_output_name = function() {
      return inputs.treatmentfile.basename.split('.').slice(0,-1).join('.')+'_p'+inputs.pvalue.split('-').slice(1)+'_kd-'+inputs.keep_dup;
   };
  - var var_output_folder = function() {
      return inputs.treatmentfile.basename.split('.').slice(0,-1).join('.')+'-p'+inputs.pvalue.split('-').slice(1)+'_kd-'+inputs.keep_dup;
   };
  - var var_output_prefix = function() {
      if (inputs.outputfolder == "") {
        if (inputs.outputname == "") {
          return inputs.treatmentfile.basename.split('.').slice(0,-1).join('.')+'-p'+inputs.pvalue.split('-').slice(1)+'_kd-'+inputs.keep_dup+'/'+inputs.treatmentfile.basename.split('.').slice(0,-1).join('.')+'_p'+inputs.pvalue.split('-').slice(1)+'_kd-'+inputs.keep_dup;
        }
        else {
          return inputs.treatmentfile.basename.split('.').slice(0,-1).join('.')+'-p'+inputs.pvalue.split('-').slice(1)+'_kd-'+inputs.keep_dup+'/'+inputs.outputname;
        }
      }
      else {
        if (inputs.outputname == "") {
          return inputs.outputfolder+'/'+inputs.treatmentfile.basename.split('.').slice(0,-1).join('.')+'_p'+inputs.pvalue.split('-').slice(1)+'_kd-'+inputs.keep_dup;
        }
        else {
          return inputs.outputfolder+'/'+inputs.outputname;
        }
      }
    };


inputs:
  treatmentfile:
    type: File
    label: "BAM file"
    inputBinding:
      prefix: -t
      position: 1

  pvalue:
    type: string?
    label: "Pvalue cutoff for peak detection"
    default: '1e-9'
    inputBinding:
      prefix: -p

  keep_dup:
    type: string?
    label: "Control towards duplicate tags"
    default: auto
    inputBinding:
      prefix: --keep-dup=
      separate: false

  wiggle:
    type: boolean?
    label: "Generate wiggle file"
    default: true
    inputBinding:
      prefix: -w

  single_profile:
    label: "Only one file for the whole genome is saved"
    type: boolean?
    default: true
    inputBinding:
      prefix: -S

  space:
    type: int?
    label: "Resolution for saving wiggle files"
    default: 50
    inputBinding:
      prefix: --space=
      separate: false

  outputname:
    type: string?
    label: "Output file prefix"
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
    label: "Output directory name"
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
    label: "Move files to new directory"
    inputBinding:
      position: 1000 
      shellQuote: false
      prefix: '" && mkdir -p $nameoffolder &&'
    default: true

  verifymove2:
    type: string?
    label: "Continuation of verifymove"
    inputBinding:
      position: 1001
      shellQuote: false
      separate: false
      prefix: 'mv '
      valueFrom: |
        ${
            if (self == ""){
              return inputs.treatmentfile.basename.split('.').slice(0,-1).join('.')+'_p'+inputs.pvalue.split('-').slice(1)+'_kd-'+inputs.keep_dup+'* $nameoffolder';
            } else {
              return self;
            }
        }
    default: ""


outputs:
  peaksbedfile:
    type: File
    label: "Peaks BED file"
    outputBinding:
      glob: |
        ${
          return var_output_prefix() + '_peaks.bed';
        }

  peaksxlsfile:
    type: File
    label: "Peaks XLS file"
    outputBinding:
      glob: |
        ${
          return var_output_prefix() + '_peaks.xls';
        }

  summitsfile:
    type: File
    label: "Peaks Summit BED file"
    outputBinding:
      glob: |
        ${
          return var_output_prefix() + '_summits.bed';
        }

  wigfile:
    type: File
    label: "Peaks WIG file"
    outputBinding:
      glob: |
        ${
          return var_output_prefix() + '_MACS_wiggle/treat/*_treat_afterfiting_all.wig.gz';
        }

  macsDir:
    type: Directory
    label: "Output directory"
    outputBinding:
      glob: |
        ${
          if (inputs.outputfolder == "") {
            return var_output_folder();
          } else {
            return inputs.outputfolder;
          }
        }
