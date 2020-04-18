# (S)ingle (E)nd (A)ntibody (SEQ)uencing pipeline

Chromatin Single-End analysis pipeline

The SEASEQ Pipeline is a complete analysis pipeline for CHiP 
sequencing single-end data.

The analysis pipeline includes mapping using bowtie, peak-calls 
using MACS and SICER, motif analysis using meme suite,
enhancers & super-enhancers using ROSE, bam density plots 
using BAM2GFF.


## INSTRUCTIONS

To run SEASEQ pipeline, you will need Linux, or some compatible 
container technology, CWL (Common Workflow Language), 
and about 30GB of supplemental data. 


## PROGRAMS & VERSIONS

* bowtie v. 1.2.2
* fastqc v. 0.11.5
* samtools v. 1.9
* R v. 3.6.1
* macs v. 041014
* SICER2 v. 1.0.1
* meme v. 4.11.2
* phantompeakqualtools v. 1.2.1.1
* bedtools/2.25.0
* python v. 3.7.0
* java v. 1.8.0_60
* perl v. 5.10.1
* wigToBigWig v. 4
* bedops v. 2.4.2
* igvtools v. 2.3.2
* ROSE v. 1.1.0
* BAM2GFF v. 1.1.0


## REQUIREMENTS

INPUT YML

```
reference: 
  class: Directory
  location: /path/to/genome_reference+index

fastqfile: 
  - { class: File, path: /path/to/fastqfile1 }
  - { class: File, path: /path/to/fastqfile2 }

keep_dup: all

chromsizes: 
  class: File
  path: /path/to/chromsizes_file

blacklistfile:
  class: File
  path: /path/to/blacklist_file

gtffile:
  class: File
  path: /path/to/gtf_file

motifdatabases:
  - { class: File, path: /path/to/meme_motif1 }
  - { class: File, path: /path/to/meme_motif2 }
```


## EXAMPLE

We provided example instructions for running under [Toil]
(https://toil.readthedocs.io/en/latest/) on our St. Jude HPC LSF cluster.


