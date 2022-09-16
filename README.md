<p align="center">
  <h1 align="center">
    SEAseq
  </h1>
  <p align="center">
   <a href="https://github.com/stjude/seaseq" target="_blank">
     <img alt="Status"
          src="https://img.shields.io/badge/status-active-success.svg" />
   </a>
   <a href="https://github.com/stjude/seaseq/issues" target="_blank">
     <img alt="Github Issues"
          src="https://img.shields.io/github/issues/stjude/seaseq" />
   </a>
   <a href="https://github.com/stjude/seaseq/pulls" target="_blank">
     <img alt="Pull Requests"
          src="https://img.shields.io/github/issues-pr/stjude/seaseq" />
   </a>
  </p>
</p>

**Single-End Antibody SEQuencing pipeline** (abbreviated as **SEAseq**) is a comprehensive automated pipeline for ChIP-Seq/CUT&RUN data analysis. Speaking broadly, it containerizes and joins field-standard, open-source tools for processing raw data and performing a wide array of basic analyses.

SEAseq analyses include alignment, peak calling, motif analysis, read coverage profiling, clustered peak (e.g. super-enhancer) identification, and quality assessment metrics, as well as automatic interfacing with data in [GEO]/[SRA]. The easy-to-use and flexibility of SEAseq makes it a reliable and efficient resource for ensuring high quality ChIP-Seq analysis, especially 
in research environments lacking computational infrastructure or expertise.  

[GEO]: https://www.ncbi.nlm.nih.gov/geo
[SRA]: https://www.ncbi.nlm.nih.gov/sra

<p align="center">
  <a href="https://github.com/stjude/seaseq/tree/master/docs#readme"><strong>Explore the documentation »</strong></a>
  <br />
  <a href="https://doi.org/10.1186/s12859-022-04588-z" target="_blank"><strong>Read the paper »</strong></a>
  <br />
  <br />
  <a href="https://github.com/stjude/seaseq/issues/new?assignees=&labels=&template=feature_request.md&title=Descriptive%20Title&labels=enhancement">Request Feature</a>
   |
  <a href="https://github.com/stjude/seaseq/issues/new?assignees=&labels=&template=bug_report.md&title=Descriptive%20Title&labels=bug">Report Bug</a>
  <br />
  ⭐ Consider starring the repo! ⭐
  <br />
</p>

## What's new in Version [3.0](https://github.com/stjude/seaseq/releases/tag/3.0)

- **PEAseq pipeline.** 
  > **PEAseq (Paired-End Antibody Sequencing Pipeline)** performs all analysis provided in SEAseq and also in a paired-end aware manner; results from SEAseq will be stored under `/single-end_mode`.

- A new color-rank scheme for the Quality Metrics and Evaluation Report HTML.
<p align="center"><img src="https://github.com/stjude/seaseq/blob/master/docs/images/colorscale.png"></p>


## SEAseq on Linux or HPC

SEAseq pipeline requires the [Cromwell] jar runner, [Docker] or [Singularity], [Java (v1.8.0)] and about 30GB of supplemental data.

View [`/test`](https://github.com/stjude/seaseq/tree/master/test) folder for example usage.

#### NOTE : HPC platforms using Singularity will require a configuration file to properly execute cromwell. Please consult [hpc-configurations](docs/hpc-configurations#readme) for more details.

[Docker]: https://www.docker.com
[Singularity]: https://docs.sylabs.io/guides/3.0/user-guide/quick_start.html
[Java (v1.8.0)]: https://java.com/en/download/help/download_options.xml
[Cromwell]: https://github.com/broadinstitute/cromwell/releases

**St. Jude Users** please consult [SEAseq on St. Jude HPC](https://wiki.stjude.org/display/compbio/sjcb+SEAseq)
## SEAseq on St. Jude cloud

Before you can run SEAseq on St. Jude Cloud, you must first create a workspace in
DNAnexus for the run. Refer to [the general workflow
guide](https://university.stjude.cloud/docs/genomics-platform/analyzing-data/running-sj-workflows#getting-started) to learn
how to create a DNAnexus workspace for each workflow run.

You can navigate to the SEAseq workflow page
[here](https://platform.stjude.cloud/workflows/seaseq).


## Citation 
Adetunji, M.O., Abraham, B.J. SEAseq: a portable and cloud-based chromatin occupancy analysis suite. BMC Bioinformatics 23, 77 (2022). https://doi.org/10.1186/s12859-022-04588-z
