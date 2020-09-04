# (S)ingle (E)nd (A)ntibody (SEQ)uencing pipeline
### SEASEQ pipeline written in WDL

Chromatin Single-End analysis pipeline

The SEASEQ Pipeline is a complete analysis pipeline for CHiP 
sequencing single-end data.

The analysis pipeline includes mapping using bowtie, peak-calls 
using MACS and SICER, motif analysis using meme suite,
enhancers & super-enhancers using ROSE, bam density plots 
using BAM2GFF.

## INSTRUCTIONS

To run SEASEQ pipeline, you will need Linux,
[Cromwell](https://github.com/broadinstitute/cromwell/releases) runner, docker,
and about 30GB of supplemental data. 

## PROGRAMS & VERSIONS

Programs and versions used to build and test the pipeline.
Programs are dockerized and do not require installation.

* bowtie v. 1.2.3
* fastqc v. 0.11.5
* samtools v. 1.9
* R v. 3.6.1
* macs v. 041014
* SICER2 v. 1.0.2
* meme v. 5.1.1
* spp v. 1.16.0
* bedtools/2.25.0
* python v. 3.7.0
* java v. 1.8.0_60
* perl v. 5.10.1
* wigToBigWig v. 4
* bedops v. 2.4.2
* igvtools v. 2.3.2
* ROSE v. 1.1.0
* BAM2GFF v. 1.1.0

## USAGE

```
java -jar cromwell.jar run seaseq.wdl -i inputs.json -o options.json
```
View /test folder for example usage and further instructions.
