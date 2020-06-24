#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [ summaryfacts.pl ]
class: CommandLineTool
label: QC on the peaks file
doc: |
  perl summaryfacts.pl -b <bam file> -p <peaks bed> -bamflag <bamflagstat> -rmdupflag <rmdupflagstat> -bkflag <bklistflagstat> -fqc <fastqczipfile> -outfile <outputfilename>


hints:
  DockerRequirement:
    dockerPull: madetunj/seaseq:v0.0.1


requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var var_output_name = function() {
      if (inputs.outfile == ""){
        if (inputs.bamflag != null) { return inputs.bamflag.nameroot+'-stats.out'; }
        else if (inputs.rmdupflag != null) { return inputs.rmdupflag.nameroot+'-stats.out'; }
        else if (inputs.bkflag != null) { return inputs.bkflag.nameroot+'-stats.out'; }
      }
   };


inputs:
  bambed:
    type: File?
    label: "Peaks BAM to BED file"
    inputBinding:
      prefix: -b

  sppfile:
    type: File?
    label: "Phantom peaks run SPP file"
    inputBinding:
      prefix: -s

  countsfile:
    type: File?
    label: "Count of peaks overlap between BED and BAMtoBED file"
    inputBinding:
      prefix: -c

  peaksxls:
    type: File?
    label: "MACS Peaks xls file"
    inputBinding:
      prefix: -px

  bamflag:
    type: File?
    label: "BAM FlagStat file"
    inputBinding:
      prefix: '-bamflag'

  rmdupflag:
    type: File?
    label: "Remove Duplicates BAM FlagStat file"
    inputBinding:
      prefix: '-rmdupflag'

  bkflag:
    type: File?
    label: "Blacklist BAM FlagStat file"
    inputBinding:
      prefix: '-bkflag'

  fastqczip:
    type: File?
    label: "FastQC ZIP file"
    inputBinding:
      prefix: '-fqc'

  fastqmetrics:
    type: File?
    label: "Basic Metrics for FastQ reads"
    inputBinding:
      prefix: '-fx'

  rosedir:
    type: Directory?
    label: "ROSE output directory"
    inputBinding:
      prefix: '-rose'

  outfile:
    type: string?
    label: "Output file name"
    inputBinding:
      position: 1000
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
    label: "Rudimentary output file"
    outputBinding: 
      glob: '*stats.out'

  htmlfile:
    type: File
    label: "HTML file"
    outputBinding:
      glob: '*stats.html'

  textfile:
    type: File
    label: "TAB-delimited file"
    outputBinding:
      glob: '*stats.txt'
