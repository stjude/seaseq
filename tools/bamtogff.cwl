#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [BAM2GFF-local.sh]
class: CommandLineTool
label: BAM to GFF for MetaGenes calculation v1 on bam file for all metagenes
doc: |
  BAM2GFF_call.sh <gtf file> <feature type> <bam file> <chromsizes file> <samplename>

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var var_output_name = function() {
      if (inputs.samplename == ""){
        return inputs.bamfile.nameroot;
      }
    };

inputs:
  bamfile:
    type: File
    inputBinding:
      position: 3
    secondaryFiles: 
      - .bai

  gtffile:
    type: File
    inputBinding:
      position: 1

  chromsizes:
    type: File
    inputBinding:
      position: 4

  feature:
    type: string?
    default: "gene"
    inputBinding:
      position: 2

  samplename:
    type: string?
    inputBinding:
      position: 5
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
    default: "bamdensity_out"

  verifymove:
    type: boolean?
    inputBinding:
      position: 1000 
      shellQuote: false
      prefix: '" && mkdir -p $nameoffolder && mv matrix *png *pdf $nameoffolder' 
    default: true

outputs:
  metagenesDir:
    type: Directory
    outputBinding:
      glob: |
        ${
          return inputs.outputfolder;
        }
