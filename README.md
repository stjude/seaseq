## Overview

**Single-End Antibody SEQuencing pipeline** (abbreviated as **SEASEQ**) is a comprehensive automated pipeline for ChIP-Seq/CUT&RUN data analysis. Speaking broadly, it containerizes and joins field-standard, open-source tools for processing raw data and performing a wide array of basic analyses.  


SEASEQ analyses include alignment, peak calling, motif analysis, read coverage profiling, clustered peak (e.g. super-enhancer) identification, and quality assessment metrics, as well as automatic interfacing with data in [GEO]/[SRA]. The easy-to-use and flexibility of SEASEQ makes it a reliable and efficient resource for ensuring high quality ChIP-Seq analysis, especially in research environments lacking computational infrastructure or expertise.  


[SRA]: https://www.ncbi.nlm.nih.gov/sra
[GEO]: https://www.ncbi.nlm.nih.gov/geo/

## Inputs

|   Name                        |   Type                | Description                                                 |  Example            |
|-------------------------------|-----------------------|-------------------------------------------------------------|---------------------|
|   FASTQ files                 |   Array of files      |   One or more FASTQ files. The files should be gzipped.     |   [`*.gz`]          |
|   SRA identifiers (SRR)       |   Array of strings    |   One or more SRRs.                                         |   [`SRR12345678`]   |
|   Genome Reference            |   File                |   The genome reference.                                     |   [`*.fa`]          |
|   Genome Bowtie indexes       |   Array of files      |   The genome bowtie v1 indexes. Should be six index files.  |   [`*.ebwt`]        |
|   Gene Annotation             |   File                |   A database file of known gene positions.                  |   [`*.gtf`]         |
|   Blacklists                  |   File                |   UHS/DER/DAC or custom blacklist regions file.             |   [`*/bed`]         |
|   MEME motif databases        |   Array of files      |   One or more of the [MEME suite motif databases]           |   [`*.meme`]        |

### Input configuration

SEASEQ supports FASTQ files and SRA identifiers (SRRs) as input. A combination of both is also supported. 
Bowtie genomic indexes and region-based blacklists are optional.

[MEME suite motif databases]: https://meme-suite.org/meme/db/motifs

## Outputs

SEASEQ provides multiple outputs from the different analysis offerings. 
Outputs are grouped into subdirectories:

| Name                  | Type    | Description                                                                            |
|-----------------------|---------|----------------------------------------------------------------------------------------|
|  [BAM_files]          | Folder  | All mapping (`.bam`) files.                                                            |
|  [BAM_Density]        | Folder  | Reads Coverage profiling in promoter and genic regions matrices, plots and heatmaps.   |
|  [MOTIFS]             | Folder  | Motifs discovery and prediction results files.                                         |
|  [PEAKS]              | Folder  | Identified narrow peaks, broad peaks and linear-stitched peaks files.                  |
|  [PEAKS_Annotation]   | Folder  | Genic annotation of peaks tables and plot.                                             |
|  [PEAKS_Display]      | Folder  | Normalized signal data tracks in wiggle, tdf and bigwig formats.                       |
|  [QC]                 | Folder  | Quality statistics and metrics of FASTQs and peaks as tables and color-coded HTML.     |

[BAM_Density]: #reads-coverage-profiling
[BAM_files]: #reads-alignment-and-filtering
[MOTIFS]: #discovery-and-enrichment-of-motifs
[PEAKS]: #peaks-identification
[PEAKS_Annotation]: #peaks-annotation
[PEAKS_Display]: #peaks-display
[QC]: #qc-metrics

## Workflow Steps
1. If provided, SRRs are downloaded as FASTQs using the [SRA Toolkit].
2. Sequencing FASTQs are aligned to the reference genome using [Bowtie].
3. Mapped reads are further processed by removal of redundant reads and blacklisted regions. 
4. Read density profiling in relevant genomic regions such as promoters and gene body using [BamToGFF].
5. Normalized and unnormalized coverage files for display are generated for the [UCSC genome browser] and [IGV]. 
6. Identification of enriched regions for two binding profiles:
    * For factors that bind shorter regions, e.g. many sequence-specific transcription factors using [MACS].
    * For broad regions of enrichment, e.g. some histone modifications using [SICER].
7. Identification of stitched clusters of enriched regions and separates exceptionally large regions, e.g. super-enhancers from typical enhancers, using [ROSE].
8. Motif discovery and enrichment using tools from the [MEME Suite].
9. Custom annotation of peaks in genic regions.
10. Assessment of quality by calculating relevant metrics including those recommmended by the [ENCODE consortium]. More information is provided [here](#qc-metrics).

## Creating a workspace

Before you can run one of our workflows, you must first create a workspace in
DNAnexus for the run. Refer to [the general workflow
guide](running-sj-workflows.md#getting-started) to learn how to create a
DNAnexus workspace for each workflow run.

You can navigate to the SEASEQ workflow page
[here](https://platform.stjude.cloud/workflows/seaseq).

## Uploading Input Files

SEASEQ requires at least the genome reference sequence, gene annotation and motif database 
files to be uploaded as [input](#inputs).

Refer to [the general workflow
guide](running-sj-workflows.md#uploading-files) to learn how to upload input
files to the workspace you just created.

## Running the Workflow

Refer to [the general workflow
guide](running-sj-workflows.md#running-the-workflow) to learn how to launch
the workflow, hook up input files, adjust parameters, start a run, and
monitor run progress.

SEASEQ has special preset options in the "Launch Tool" dropdown.
These are explained below:

## Analysis of Results

Each tool in St. Jude Cloud produces a visualization that makes understanding 
results more accessible than working with spreadsheets or tab-delimited files. 
This is a recommended way to view the provided visualization files.

Refer to [the general workflow
guide](running-sj-workflows.md#custom-visualizations) to learn how to access
these visualizations.

We also include the raw output files for you to dig into if the visualization 
is not sufficient to answer your research question.

Refer to [the general workflow
guide](running-sj-workflows.md#raw-results-files) to learn how to access raw
results files.

SEASEQ results will be in the parent `/` folder unless otherwise specified. 

##  Interpreting results

Upon successful run of SEASEQ, all files are saved to the results directory as 
listed in the [Outputs Section](#outputs). Here, we will discuss some of the
different output directories in more detail.


### Reads Alignment and Filtering

Reads are stingently mapped to reference genome using [Bowtie] 
`bowtie -l readlength -p 20 -k 2 -m 2 --best -S`. 

Mapped BAMs are further processed by removal of duplicate reads using [SAMTools] 
and blacklisted regions using [bedtools].
Blacklisted regions are sources of bias in most ChIP-Seq experiments. These problematic regions are identified as regions with significant background noise or artifically high signal (UHS/DAC/DER), and are recommended to be excluded in order to assess biologically relevant and true signals of enrichment.

### Peaks Identification

We identify enriched regions based on two binding profiles; factors that bind 
shorter regions such as transcription factors using MACS, and for those that 
bind broad regions such as some histone marks using SICER. 

Parameters specified for peak calls : 
  * **Narrow Peaks** : `macs14 -p 1e-9 --keep-dup=auto --space=50 -w -S`
  * **Broad Peaks** : `sicer -rt 1 -w 200 -f 150 -egf 0.86 -g 200 -e 100`

In addition, we identify large clusters of enrichment (enhancers) and 
exceptionally large regions (super-enhancers) using the [ROSE] algorithm.
Computed stitched regions are generated as tab-delimited or BED files with 
overlapping gene information. 

### Peaks Display

Coverage graphical files are generated in different visualization formats; 
[WIG], [bigWig] and [TDF] formats, and normalized for display on 
multiple genomic browsers.

Parameter specified to generate coverage files: `macs14 --space=50 --shiftsize=200 -w -S`.

### Discovery and Enrichment of Motifs

We discover sequence patterns that are widespread and have biological 
significance using the [MEME-ChIP] and [AME] tools from the [MEME Suite].

[AME] discovers known enriched binding motifs. 
[MEME-ChIP] performs several motif analysis steps such as novel DNA-binding motif discovery 
(MEME and DREME), identify motif enrichment relative to background (CentriMo), visualize 
the arrangement of the predicted motif sites (SpaMO, FIMO).

The motifs are identified using the peak regions and 100bp window around peak summit (-50 and +50).

### Reads Coverage Profiling

Read density profiling of major genomic regions such as promoters and gene body using [BamToGFF].
We generate coverage graphs and heatmaps plots with additional custom R scripts for further customization.

### Peaks Annotation

Genic annotation of peaks including promoters, gene bodies, gene-centric windows, and proximal genes.
We designed custom scripts to provide this information.

Custom scripts are designed to generate the genic annotation of peaks at promoters, gene bodies 
and gene-centric windows. Annotated regions are collated to provide a binary or matrix overview 
of proximal genes, distribution graphs are also provided.

### QC Metrics

SEASEQ provides a vast set of quality metrics for detecting experimental issues, including ChIP-Seq 
metrics recommended by the ENCODE consortium. 

We incorporated a five-scale color-rank flag system to visually identify excellent 
(score = 2), good (score = 1), average (score = 0), below-average (score = -1) or poor (score = -2) 
results for each metric and a cross-metric summary score (between -2 and 2), using recommended 
thresholds where possible. 

The metrics are color flagged for easy visualization of overall performance in HTML format.

SEASEQ metrics calculated to infer quality are:
| Quality Metric	| Definition	|
|--|--|
| Aligned Percent | Percentage of mapped reads. |
| Base Quality	| Per-base sequence quality. |
| Estimated Fragment Width	| Average fragment size of the peak distribution.	|
| Estimated Tag Length	| Sequencing read length. |
| FRiP	| The fraction of reads within peaks regions. |
| Linear Stitched Peaks (Enhancers)	| Total number of clustered enriched regions. |
| Non-Redundant Fraction (NRF)	| Fraction of uniquely mapped reads. |
| Normalized Strand-correlation Coefficient (NSC)	| To determine signal-to-noise ratio using strand cross-correlation. The ratio of the maximum cross-correlation value divided by the background cross-correction. |
| Sequence Diversity	| Sequence overrepresentation.  If reads/sequences are overrepresented in the library. |
| PCR Bottleneck Coefficient (PBC)	| It is a measure of library complexity determined by the fraction of genomic locations with exactly one unique read versus those covered by at least one unique reads.	|
| Peaks	| Total number of enriched regions. |
| Raw reads	| Total number of sequencing reads. |
| Relative Strand-correlation Coefficient (RSC)	| A strand cross-correlation ratio between the fragment-length cross-correlation and the read-length peak. |
| SE-like enriched regions (Super Enhancers)	| Total number of SE-like clustered enriched regions. |
| Overall Quality	| Cross-metric average score. |

## Frequently asked questions

None yet!

If you have any questions not covered here, feel free to reach out
on [our contact
form](https://hospital.stjude.org/apps/forms/fb/st-jude-cloud-contact/).

# References

None yet!

# Similar Topics

[Running our Workflows](../analyzing-data/running-sj-workflows.md)

[Working with our Data Overview](../managing-data/working-with-our-data.md)

[Upload/Download Data (local)](../managing-data/upload-local.md)



[SAMTools]: https://doi.org/10.1093/bioinformatics/btp352
[bedtools]: https://doi.org/10.1093/bioinformatics/btq033
[SRA Toolkit]: http://www.ncbi.nlm.nih.gov/books/NBK158900/
[Bowtie]: https://doi.org/10.1186/gb-2009-10-3-r25
[BamToGFF]: https://github.com/stjude/BAM2GFF
[UCSC genome browser]: https://doi.org/10.1093/bib/bbs038
[IGV]: https://doi.org/10.1093/bib/bbs017
[MACS]: https://doi.org/10.1186/gb-2008-9-9-r137
[SICER]: https://doi.org/10.1093/bioinformatics/btp340
[ROSE]: http://younglab.wi.mit.edu/super_enhancer_code.html
[MEME Suite]: https://doi.org/10.1093/nar/gkv416
[ENCODE consortium]: https://doi.org/10.1101/gr.136184.111
[WIG]: https://genome.ucsc.edu/goldenPath/help/wiggle.html 
[bigWig]: https://genome.ucsc.edu/goldenpath/help/bigWig.html
[TDF]: https://software.broadinstitute.org/software/igv/TDF 
[AME]: https://meme-suite.org/meme/tools/ame
[MEME-ChIP]: https://meme-suite.org/meme/tools/meme-chip
[seaseq]: https://github.com/stjude/seaseq
