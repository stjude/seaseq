#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: bowtie
class: CommandLineTool

#INITIAL SYNTAX
label: Bowtie on ChipSeq reads
doc: |
  bsub -K -R \"rusage[mem=10000]\" -n 20 -q standard -R \"select[rhel7]\" bowtie -l \$readLength -p 20 -k 2 -m 2 --best -S /datasets/public/genomes/Homo_sapiens/hg19/TOPHAT/hg19 ${files[$file]} \> $file.sam" >> $outfile


requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

  expressionLib:
  - var var_output_name = function() {
      if (inputs.output_prefix != null) {
        return inputs.output_prefix+'.sam';
      } else {
        return inputs.fastqfile.basename.split('.fastq')[0]+'.sam';
      }
   };
  - var var_readLength = function() {
      if (inputs.readLengthFile != null) {
        return inputs.readLengthFile.nameroot;
      }
   };

inputs:
  output_prefix:
    type: string?

  processors:
    type: int?
    default: 20
    inputBinding:
      prefix: -p
      position: 1

  good_alignments:
    type: int?
    default: 2
    inputBinding:
      prefix: -k
      position: 1

  limit_alignments:
    type: int?
    default: 2
    inputBinding:
      prefix: -m
      position: 1

  readLength:
    type: int
    inputBinding:
      prefix: -l
      position: 2
      valueFrom: |
        ${
            if (self == 0){
              return var_readLength();
            } else {
              return self;
            }
         }
    default: 0

  readLengthFile:
    type: File?

  best_alignments:
    type: boolean?
    default: true
    inputBinding:
      prefix: --best
      position: 3

  reference:
    type: Directory
    inputBinding:
      position: 4
      valueFrom: |
        ${
            for (var i = 0; i < self.listing.length; i++) {
                if (self.listing[i].path.split('.').slice(-3).join('.') == 'rev.1.ebwt') {
                  return self.listing[i].path.split('.').slice(0,-3).join('.');
                }
            }
            return null;
        }
    doc: |
      Folder with Bowtie indices

  fastqfile:
    type: File
    inputBinding:
      position: 5

  samoutput:
    type: boolean?
    default: true
    inputBinding:
      prefix: -S
      position: 6

  samfile:
    type: string
    inputBinding:
      prefix: '>'
      position: 7
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
  samfile:
    type: File
    outputBinding:
      glob: '*sam'
