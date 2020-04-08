#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [ mkdir, -p ]
class: CommandLineTool 

label: move files

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var var_output_name = function() {
      if (inputs.readqc_zip != null) { return inputs.readqc_zip.basename.split("_fastqc").slice(0,-1); }
   };


inputs:
  makedirname:
    type: string?
    inputBinding:
      position: 1
      shellQuote: false
      valueFrom: |
        ${
            if (self == ""){
              return var_output_name();
            } else {
              return self;
            }
        }
    default: ""   
  
  copyoutput:
    type: string?
    inputBinding:
      position: 3
      shellQuote: false
      prefix: '&& cp -rf'
    default: ""
    
  sam_sort:
    type: File?
    inputBinding:
      position: 10

  fastq_metrics:
    type: File?
    inputBinding:
      position: 10
      
  rmdup_bam:
    type: File?
    inputBinding:
      position: 10
      
  bklist_bam: 
    type: File?
    inputBinding:
      position: 10
      
  bamqc_html:
    type: File?
    inputBinding:
      position: 10
      
  bamqc_zip:
    type: File?
    inputBinding:
      position: 10
      
  readqc_zip:
    type: File
    inputBinding:
      position: 10
      
  readqc_html:
    type: File?
    inputBinding:
      position: 10
      
  macsDir:
    type: Directory?
    inputBinding:
      position: 10
      
  allmacsDir:
    type: Directory?
    inputBinding:
      position: 10
      
  nmmacsDir:
    type: Directory?
    inputBinding:
      position: 10
      
  rpmwig:
    type: File?
    inputBinding:
      position: 10
      
  outBW:
    type: File?
    inputBinding:
      position: 10
      
  outtdf:
    type: File?
    inputBinding:
      position: 10
      
  allrpmwig:
    type: File?
    inputBinding:
      position: 10
      
  alloutBW:
    type: File?
    inputBinding:
      position: 10
      
  allouttdf:
    type: File?
    inputBinding:
      position: 10
      
  nmrpmwig:
    type: File?
    inputBinding:
      position: 10
      
  nmoutBW:
    type: File?
    inputBinding:
      position: 10
      
  nmouttdf:
    type: File?
    inputBinding:
      position: 10
      
  bedfasta:
    type: File?
    inputBinding:
      position: 10
      
  flankbed:
    type: File?
    inputBinding:
      position: 10
      
  memechipdir:
    type: Directory?
    inputBinding:
      position: 10
      
  summitmemechipdir:
    type: Directory?
    inputBinding:
      position: 10
      
  amedir:
    type: Directory?
    inputBinding:
      position: 10
      
  summitamedir:
    type: Directory?
    inputBinding:
      position: 10
      
  metagenesDir:
    type: Directory?
    inputBinding:
      position: 10
      
  sicerDir:
    type: Directory?
    inputBinding:
      position: 10
      
  roseoutput:
    type: Directory?
    inputBinding:
      position: 10
      
  statsfile:
    type: File?
    inputBinding:
      position: 10
      
  htmlfile:
    type: File?
    inputBinding:
      position: 10
      
  textfile:
    type: File?
    inputBinding:
      position: 10
      
  outputdir:
    type: string?
    inputBinding:
      position: 1000
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
  finalDir:
    type: Directory
    outputBinding:
      glob: |
        ${
          if (inputs.outputdir == "") {
            return var_output_name();
          } else {
            return inputs.outputdir;
          }
        }
