# SEAseq Quality Metrics Expatiated Terms and Definitions

To view the complete list of metrics SEAseq offersMost metrics are adopted 
from [ENCODE]([Landt et al, Genome Res.2012]).

Additional information on applicable ChIP-seq metrics are provided 
[here](https://genome.ucsc.edu/ENCODE/qualityMetrics.html).

[ENCODE]: https://www.encodeproject.org/data-standards/terms
[Landt et al, Genome Res.2012]: https://doi.org/10.1101/gr.136184.111


## Fraction of Reads in Peaks (FRiP)

Fraction of all mapped reads that fall into the called peak 
regions, i.e. usable reads in significantly enriched peaks 
divided by all usable reads. In general, FRiP scores 
correlate positively with the number of regions.

## Non-Redundant Fraction (NRF)

Number of distinct uniquely mapping reads 
(i.e. after removing duplicates) / Total number of reads.

## Normalized Strand-correlation Coefficient (NSC)

The NSC is the ratio of the maximal cross-correlation value 
(which occurs at strand shift equal to fragment length) 
divided by the background cross-correlation (minimum 
cross-correlation value over all possible strand shifts).

## Relative Strand-correlation Coefficient (RSC)

The RSC is the ratio of the fragment-length cross-correlation 
value minus the background cross-correlation value, divided by 
the phantom-peak cross-correlation value minus the background 
cross-correlation value. The minimum possible value is 0 
(no signal), highly enriched experiments have values greater 
than 1, and values much less than 1 may indicate low quality.

## PCR Bottleneck Coefficient (PBC)

A measure of library complexity, i.e. 
how skewed the distribution of read 
counts per location is towards 1 read per location.

PBC = N1/Nd

(where N1= number of genomic locations to which EXACTLY 
one unique mapping read maps, and Nd = the number of 
genomic locations to which AT LEAST one unique mapping 
read maps, i.e. the number of non-redundant, unique mapping reads).


### Disclaimer

All definitions were copied from the resources cited below:
1. [ENCODE Quality Metrics](https://genome.ucsc.edu/ENCODE/qualityMetrics.html).
1. [ENCODE Data Standards](https://www.encodeproject.org/data-standards/terms).
1. [Landt et al, Genome Res.2012](https://doi.org/10.1101/gr.136184.111).
