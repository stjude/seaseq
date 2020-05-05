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
    label: "Make directory"
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
    label: "copy all files"
    inputBinding:
      position: 3
      shellQuote: false
      prefix: '&& cp -rf'
    default: ""
    
  sam_sort:
    type: File?
    label: "BAM file"
    inputBinding:
      position: 10

  fastq_metrics:
    type: File?
    label: "Basic FastQ Metrics file"
    inputBinding:
      position: 10
      
  rmdup_bam:
    type: File?
    label: "Remove duplicates BAM file"
    inputBinding:
      position: 10
      
  bklist_bam: 
    type: File?
    label: "Blacklist BAM file"
    inputBinding:
      position: 10
      
  bamqc_html:
    type: File?
    label: "BAM fastQC HTML file"
    inputBinding:
      position: 10
      
  bamqc_zip:
    type: File?
    label: "BAM fastQC ZIP file"
    inputBinding:
      position: 10
      
  readqc_zip:
    type: File
    label: "FastQ fastQC ZIP file"
    inputBinding:
      position: 10
      
  readqc_html:
    type: File?
    label: "FastQ fastQC HTML file"
    inputBinding:
      position: 10
      
  macsDir:
    type: Directory?
    label: "MACS Output directory"
    inputBinding:
      position: 10
      
  allmacsDir:
    type: Directory?
    label: "MACS ALL Output directory"
    inputBinding:
      position: 10
      
  nmmacsDir:
    type: Directory?
    label: "MACS nomodel output directory"
    inputBinding:
      position: 10
      
  rpmwig:
    type: File?
    label: "MACS wig file"
    inputBinding:
      position: 10
      
  outBW:
    type: File?
    label: "MACS BigWig file"
    inputBinding:
      position: 10
      
  outtdf:
    type: File?
    label: "MACS TDF file"
    inputBinding:
      position: 10
      
  allrpmwig:
    type: File?
    label: "MACS ALL wig file"
    inputBinding:
      position: 10
      
  alloutBW:
    type: File?
    label: "MACS ALL BigWig file"
    inputBinding:
      position: 10
      
  allouttdf:
    type: File?
    label: "MACS ALL TDF file"
    inputBinding:
      position: 10
      
  nmrpmwig:
    type: File?
    label: "MACS NM wig file"
    inputBinding:
      position: 10
      
  nmoutBW:
    type: File?
    label: "MACS NM BigWig file"
    inputBinding:
      position: 10
      
  nmouttdf:
    type: File?
    label: "MACS NM TDF file"
    inputBinding:
      position: 10
      
  bedfasta:
    type: File?
    label: "BED to FASTA file"
    inputBinding:
      position: 10
      
  flankbed:
    type: File?
    label: "Flank BED file"
    inputBinding:
      position: 10
      
  memechipdir:
    type: Directory?
    label: "MEME-CHIP directory"
    inputBinding:
      position: 10
      
  summitmemechipdir:
    type: Directory?
    label: "Summit MEME-CHIP directory"
    inputBinding:
      position: 10
      
  amedir:
    type: Directory?
    label: "AME directory"
    inputBinding:
      position: 10
      
  summitamedir:
    type: Directory?
    label: "Summit AME directory"
    inputBinding:
      position: 10
      
  metagenesDir:
    type: Directory?
    label: "BAM2GFF directory"
    inputBinding:
      position: 10
      
  sicerDir:
    type: Directory?
    label: "SICER directory"
    inputBinding:
      position: 10
      
  roseoutput:
    type: Directory?
    label: "ROSE directory"
    inputBinding:
      position: 10
      
  statsfile:
    type: File?
    label: "Overall statistics file"
    inputBinding:
      position: 10
      
  htmlfile:
    type: File?
    label: "Overall statistics HTML file"
    inputBinding:
      position: 10
      
  textfile:
    type: File?
    label: "Overal statistics TAB file"
    inputBinding:
      position: 10
      
  outputdir:
    type: string?
    label: "Output directory name"
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
    label: "Output directory name"
    outputBinding:
      glob: |
        ${
          if (inputs.outputdir == "") {
            return var_output_name();
          } else {
            return inputs.outputdir;
          }
        }
