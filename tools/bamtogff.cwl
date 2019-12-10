#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [ BamToGFFMetaGenes.pl ]
class: CommandLineTool
label: BAM to GFF for MetaGenes calculation v1 on bam file for all metagenes
doc: perl BamToGFFMetaGenes.pl -g <gff file|gtf file> -b <bam file> [-o <outfile name>]

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
      prefix: -b

  gtffile:
    type: File
    inputBinding:
      prefix: -g

  chromsizes: 
    type: File
    inputBinding:
      prefix: -c

  feature:
    type: string?
    default: "gene"
    inputBinding:
      prefix: -f

  outfile:
    type: string?
    inputBinding:
      prefix: -o
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
      prefix: '&& mkdir -p bamdensity_out && mv *txt *png *pdf bamdensity_out'
    default: ""


outputs:
  metagenesDir:
    type: Directory
    outputBinding:
      glob: "bamdensity_out"

  promoters:
    type: File
    outputBinding: 
      glob: '*-promoters.pdf'

  genebody:
    type: File
    outputBinding:
      glob: '*-entiregene.pdf'

  promotersheatmap:
    type: File
    outputBinding: 
      glob: '*heatmap.promoters.png'

  genebodyheatmap:
    type: File
    outputBinding:
      glob: '*heatmap.entiregene.png'

  promotersheatmappdf:
    type: File
    outputBinding:
      glob: '*heatmap.promoters.pdf'

  genebodyheatmappdf:
    type: File
    outputBinding:
      glob: '*heatmap.entiregene.pdf'

  liqpromoters:
    type: File
    outputBinding:
      glob: '*-promoters.txt'

  liqgenes:
    type: File
    outputBinding:
      glob: '*-genebody.txt'

  liqdownstream:
    type: File
    outputBinding:
      glob: '*-downstream.txt'

  liqupstream:
    type: File
    outputBinding:
      glob: '*-upstream.txt'


