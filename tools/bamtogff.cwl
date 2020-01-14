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
      if (inputs.outfile == ""){
        if (inputs.bamfile != null) { return inputs.bamfile.nameroot; }
      }
    };

inputs:
  bamfile:
    type: File
    secondaryFiles: $(self.basename+".bai")
    inputBinding:
      position: 3

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

  savedDir:
    type: string?
    inputBinding:
      position: 1000
      shellQuote: false
      prefix: '&& mkdir -p bamdensity_out && mv matrix *png *pdf bamdensity_out'
    default: ""


outputs:
  metagenesDir:
    type: Directory
    outputBinding:
      glob: "bamdensity_out"
