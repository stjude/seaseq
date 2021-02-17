## Overview

**Single-End Analysis SEQuencing pipeline** (abbreviated as **SEASEQ**) is a comprehensive automated computational pipeline for all ChIP-Seq/CUT&RUN data analysis.

SEASEQ performs extensive analyses from the raw output of the experiment, including alignment, peak calling, motif analysis, read coverage profiling, clustered peak (e.g. super-enhancer) identification, and quality assessment metrics, as well as automatic interfacing with data in [GEO]/[SRA]. The easy-to-use and versatile format of SEASEQ makes it a reliable and efficient resource for ensuring high quality ChIP-Seq analysis.


[SRA]: https://www.ncbi.nlm.nih.gov/sra
[GEO]: https://www.ncbi.nlm.nih.gov/geo/

## Inputs

|   Name                        |   Type                | Description                                                 |  Example            |
|-------------------------------|-----------------------|-------------------------------------------------------------|---------------------|
|   FastQ files                 |   Array of files      |   One or more fastq files. The files should be gzipped.     |   [`*.gz`]          |
|   SRA identifier (SRR)        |   Array of strings    |   One of more SRRs.                                         |   [`SRR12345678`]   |
|   Genome Reference            |   File                |   The genome reference.                                     |   [`*.fa`]          |
|   Genome Bowtie indexes       |   Array of files      |   The genome bowtie v1 indexes. Should be six index files.  |   [`*.ebwt`]        |
|   Gene Annotation             |   File                |   The gene annotation.                                      |   [`*.gtf`]         |
|   Blacklists                  |   File                |   UHS/DER/DAC or custom blacklists regions file.            |   [`*/bed`]         |
|   MEME motif databases        |   Array of files      |   One or more of the [MEME suite motif databases]           |   [`*.meme`]        |

### Input configuration

SEASEQ requires either or both FastQ files or SRA identifiers (SRRs). Genome Bowtie indexes, Blacklists are optional

[MEME suite motif databases]: https://meme-suite.org/meme/db/motifs

## Output

SEASEQ provides multiple outputs from the different analyses offers. Outputs are grouped into folders

| Name                  | Type    | Description                                                                            |
|-----------------------|---------|----------------------------------------------------------------------------------------|
|  BAM_Density          | Folder  | Reads Coverage profiling in promoter and genic regions matrices, plots and heatmaps.   |
|  BAM_files            | Folder  | All mapping (.bam) files.                                                              |
|  MOTIFS               | Folder  | Motifs discovery and prediction results files.                                         |
|  PEAKS                | Folder  | Identified narrow peaks, broad peaks and linear-stitched peaks files.                  |
|  PEAKS_Annotation     | Folder  | Genic annotation of peaks tables and plot.                                             |
|  PEAKS_Display        | Folder  | Normalized signal data tracks in wiggle, tdf and bigwig formats.                       |
|  [QC]                 | Folder  | Quality statistics and metrics of fastQs and peaks as tables and color-coded HTML.     |

[QC]: #qc-metrics

## Creating a workspace

Before you can run one of our workflows, you must first create a workspace in
DNAnexus for the run. Refer to [the general workflow
guide](running-sj-workflows.md#getting-started) to learn how to create a
DNAnexus workspace for each workflow run.

You can navigate to the SEASEQ workflow page
[here](https://platform.stjude.cloud/workflows/seaseq).

## Uploading Input Files

SEASEQ requires at least the genome reference, gene annotation and motif database 
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

### Selecting parameters

The following are the parameters that can be set, a short description of each 
parameter and an example value. How to customize parameters is covered in [the 
general workflow guide](running-sj-workflows.md#selecting-parameters). If you 
have questions, please [contact us](https://stjude.cloud/contact)

| Parameter Name        | Description                                            |  Example           |
|-----------------------|--------------------------------------------------------|--------------------|
|   sra_id              |   One or more SRA identifiers (SRR)                    |   [`SRR12345678`]  |
|   gtf_feature         |   gene annotation feature (either gene or transcript)  |   [`gene`]         |

## Analysis of Results

Refer to [the general workflow
guide](running-sj-workflows.md#custom-visualizations) to learn how to access
these visualizations.

We include the output files for you to dig into to answer your research question.

Refer to [the general workflow
guide](running-sj-workflows.md#raw-results-files) to learn how to access raw
results files.

SEASEQ results will be in the parent `/` folder unless otherwise specified. 

##  Interpreting results

Upon successful run of SEASEQ, all files are saved to the results directory as specified in [outputs](#outputs)

### QC Metrics

SEASEQ provides a vast set of quality metrics for detecting experimental issues, including ChIP-Seq 
metrics recommended by the ENCODE consortium. 

SEASEQ QC metrics include the percentage of reads mapped, uniquely mapped reads (NRF), read 
quality, sequence diversity, fraction of reads in peaks (FRiP), strand correlation scores 
(NSC, RSC), library complexity (PCR bottleneck) and many other important metrics used to 
infer quality. 

In addition, we devised a five-scale color-rank flag system to visually identify excellent 
(score = 2), good (score = 1), average (score = 0), below-average (score = -1) or poor (score = -2) 
results for each metric and a cross-metric summary score (between -2 and 2), using recommended 
thresholds where possible. 
The results are exported in a tab-delimited file and color flagged in an html format.

## Frequently asked questions

None yet! If you have any questions not covered here, feel free to reach out
on [our contact
form](https://hospital.stjude.org/apps/forms/fb/st-jude-cloud-contact/).

# References

None yet!


## Similar Topics

[Running our Workflows](../analyzing-data/running-sj-workflows.md)

[Working with our Data Overview](../managing-data/working-with-our-data.md)

[Upload/Download Data (local)](../managing-data/upload-local.md)


[seaseq]: https://github.com/stjude/seaseq

