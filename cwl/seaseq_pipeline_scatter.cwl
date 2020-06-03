#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow
doc: |
  - seaseq_pipeline.cwl is a complete analysis pipeline for CHiP 
      sequencing single-end data.
  - analysis pipeline includes mapping using bowtie, peak-calls 
      using MACS1 and SICER, motif analysis using meme suite, 
      enhancers & super-enhancers using ROSE, bam density plots 
      using BAM2GFF.

requirements:
  - class: SubworkflowFeatureRequirement
  - class: ScatterFeatureRequirement

inputs:
# main files & directorys
  reference: Directory
  gtffile: File
  fastqfile: File[]
  chromsizes: File
  blacklistfile: File
  motifdatabases: File[]

# optional parameters
  #BOWTIE
  best_alignments: boolean?
  good_alignments: int?
  limit_alignments: int?
  processors: int?
  
  # MACS
  nomodel: boolean?
  wiggle: boolean?
  single_profile: boolean?
  shiftsize: int?
  space: int?
  pvalue: string?
  keep_dup: string?
  flank: int?
  
  # SICER & ROSE
  species: string?
  
  # SICER
  redundancy: int?
  window: int?
  fragment_size: int?
  genome_fraction: double?
  gapsize: int?
  evalue: int?

  # ROSE
  feature: string?

outputs:
  finalDir:
    type: Directory[]
    label: "output directory of the analysis result for each fastqfile"
    outputSource: MoveFiles/finalDir

steps:
  BasicMetrics:
    requirements:
      ResourceRequirement:
        ramMax: 20000
        coresMin: 1
    in: 
      fastqfile: fastqfile
    out: [metrics_out]
    run: basicfastqstats.cwl
    scatter: fastqfile

  TagLen:
    in: 
      datafile: BasicMetrics/metrics_out
    out: [tagLength]
    run: taglength.cwl
    scatter: datafile
   
  ReadQC:
    in:
      infile: fastqfile
    out: [htmlfile, zipfile]
    run: fastqc.cwl
    scatter: infile

  Bowtie:
    requirements:
      ResourceRequirement:
        ramMax: 10000
        coresMin: 20
    run: bowtie.cwl
    in:
      readLengthFile: TagLen/tagLength
      best_alignments: best_alignments
      good_alignments: good_alignments
      fastqfile: fastqfile
      limit_alignments: limit_alignments
      processors: processors
      reference: reference
    out: [samfile]
    scatter: [readLengthFile, fastqfile]
    scatterMethod: dotproduct

  SamView:
    in:
      infile: Bowtie/samfile
    out: [outfile]
    run: samtools-view.cwl
    scatter: infile

  BamQC:
    in:
      infile: SamView/outfile
    out: [htmlfile, zipfile]
    run: fastqc.cwl
    scatter: infile

  SamSort:
    in:
      infile: SamView/outfile
    out: [outfile]
    run: samtools-sort.cwl
    scatter: infile

  BkList:
    in:
      infile: SamSort/outfile
      blacklistfile: blacklistfile
    out: [outfile]
    run: blacklist.cwl
    scatter: infile

  BkIndex:
    in:
      infile: BkList/outfile
    out: [outfile]
    run: samtools-index.cwl
    scatter: infile

  SamRMDup:
    in:
      infile: BkList/outfile
    out: [outfile]
    run: samtools-mkdupr.cwl
    scatter: infile

  SamIndex:
    in:
      infile: SamRMDup/outfile
    out: [outfile]
    run: samtools-index.cwl
    scatter: infile

  STATbam:
    in:
      infile: SamView/outfile
    out: [outfile]
    run: samtools-flagstat.cwl
    scatter: infile

  STATrmdup:
    in:
      infile: SamRMDup/outfile
    out: [outfile]
    run: samtools-flagstat.cwl
    scatter: infile

  STATbk:
    in:
      infile: BkList/outfile
    out: [outfile]
    run: samtools-flagstat.cwl
    scatter: infile

# PEAK CALLING & VISUALS
  MACS-Auto:
    requirements:
      ResourceRequirement:
        ramMax: 10000
        coresMin: 1
    in:
      treatmentfile: BkIndex/outfile
      space: space
      pvalue: pvalue
      wiggle: wiggle
      single_profile: single_profile
    out: [ peaksbedfile, peaksxlsfile, summitsfile, wigfile, macsDir ]
    run: macs1call.cwl
    scatter: treatmentfile

  WIG-Auto:
    in:
      wigfile: MACS-Auto/wigfile
      peaksxls: MACS-Auto/peaksxlsfile
      chromsizes: chromsizes
    out: [ rpmwig, outBW, outtdf ]
    run: visualization.cwl
    scatter: [wigfile, peaksxls]
    scatterMethod: dotproduct

  MACS-All:
    requirements:
      ResourceRequirement:
        ramMax: 10000
        coresMin: 1
    in:
      treatmentfile: BkIndex/outfile
      keep_dup: keep_dup
      space: space
      pvalue: pvalue
      wiggle: wiggle
      single_profile: single_profile
    out: [ peaksbedfile, peaksxlsfile, summitsfile, wigfile, macsDir ]
    run: macs1call.cwl
    scatter: treatmentfile

  WIG-All:
    in:
      wigfile: MACS-All/wigfile
      peaksxls: MACS-All/peaksxlsfile
      chromsizes: chromsizes
    out: [ rpmwig, outBW, outtdf ]
    run: visualization.cwl
    scatter: [wigfile, peaksxls]
    scatterMethod: dotproduct

  MACS-NM:
    requirements:
      ResourceRequirement:
        ramMax: 10000
        coresMin: 1
    in:
      treatmentfile: BkIndex/outfile
      space: space
      wiggle: wiggle
      single_profile: single_profile
    out: [ peaksbedfile, peaksxlsfile, summitsfile, wigfile, macsDir ]
    run: macs1nm.cwl
    scatter: treatmentfile

  WIG-NM:
    in:
      wigfile: MACS-NM/wigfile
      peaksxls: MACS-NM/peaksxlsfile
      chromsizes: chromsizes
    out: [ rpmwig, outBW, outtdf ]
    run: visualization.cwl
    scatter: [wigfile, peaksxls]
    scatterMethod: dotproduct

# MOTIF analysis
  MOTIFS:
    in:
      reference: reference
      bedfile: MACS-Auto/peaksbedfile
      motifdatabases: motifdatabases
    out: [memechipdir, amedir, bedfasta]
    run: motifs.cwl
    scatter: bedfile

#SUMMIT-MOTIF analysis
  FlankBED:
    in:
      bedfile: MACS-Auto/summitsfile
      flank: flank
    out: [outfile]
    run: flankbed.cwl
    scatter: bedfile
  
  SummitMOTIFS:
    in:
      reference: reference
      bedfile: FlankBED/outfile
      motifdatabases: motifdatabases
    out: [memechipdir, amedir, bedfasta]
    run: motifs.cwl
    scatter: bedfile
    
# METAGENE analysis
  MetaGene:
    in:
      bamfile: SamIndex/outfile
      gtffile: gtffile
      chromsizes: chromsizes
    out:  [ metagenesDir ] 
    run: bamtogff-scatter.cwl
    scatter: bamfile

# SICER broad peaks caller
  B2Bed:
    in:
      infile: SamIndex/outfile
    out: [ outfile ]
    run: bamtobed.cwl
    scatter: infile

  SICER:
    requirements:
      ResourceRequirement:
        ramMax: 10000
        coresMin: 1
    in:
      species: species
      redundancy: redundancy
      window: window
      fragment_size: fragment_size
      genome_fraction: genome_fraction
      gapsize: gapsize
      evalue: evalue
      treatmentbedfile: B2Bed/outfile
    run: sicer.cwl
    out: [ sicerDir ]
    scatter: treatmentbedfile

# ROSE enhancer caller
  ROSE:
    requirements:
      ResourceRequirement:
        ramMax: 20000
        coresMin: 1
    in:
      species: species
      feature: feature
      gtffile: gtffile
      bamfile: BkIndex/outfile
      fileA: MACS-All/peaksbedfile
      fileB: MACS-Auto/peaksbedfile
    out: [ RoseDir ]
    run: roseNC-scatter.cwl
    scatter: [bamfile, fileA, fileB]
    scatterMethod: dotproduct

# Quality Control & Statistics
  Bklist2Bed:
    in:
      infile: BkIndex/outfile
    out: [ outfile ]
    run: bamtobed.cwl
    scatter: infile

  SortBed:
    requirements:
      ResourceRequirement:
        ramMax: 10000
        coresMin: 1
    in:
      infile: Bklist2Bed/outfile
    out: [outfile]
    run: sortbed.cwl
    scatter: infile

  runSPP:
    requirements:
      ResourceRequirement:
        ramMax: 10000
        coresMin: 1
    in:
      infile: BkIndex/outfile
    out: [spp_out]
    run: runSPP.cwl
    scatter: infile

  CountIntersectBed:
    requirements:
      ResourceRequirement:
        ramMax: 10000
        coresMin: 1
    in:
      peaksbed: MACS-Auto/peaksbedfile
      bamtobed: SortBed/outfile
    out: [outfile]
    run: intersectbed.cwl
    scatter: [peaksbed, bamtobed]
    scatterMethod: dotproduct

  PeaksQC:
    requirements:
      ResourceRequirement:
        ramMax: 10000
        coresMin: 1
    in:
      fastqmetrics: BasicMetrics/metrics_out
      fastqczip: ReadQC/zipfile
      sppfile: runSPP/spp_out
      bambed: Bklist2Bed/outfile
      countsfile: CountIntersectBed/outfile
      peaksxls: MACS-Auto/peaksxlsfile
      bamflag: STATbam/outfile
      rmdupflag: STATrmdup/outfile
      bkflag: STATbk/outfile
      rosedir: ROSE/RoseDir
    out: [ statsfile, htmlfile, textfile ]
    run: summarystats.cwl
    scatter: [fastqmetrics, fastqczip, sppfile, bambed, countsfile, peaksxls, bamflag, rmdupflag, bkflag, rosedir]
    scatterMethod: dotproduct

  MoveFiles:
    requirements:
      ResourceRequirement:
        ramMax: 10000
        coresMin: 1
    in:
      sam_sort: SamSort/outfile
      fastq_metrics: BasicMetrics/metrics_out
      rmdup_bam: SamIndex/outfile
      bklist_bam: BkIndex/outfile
      bamqc_html: BamQC/htmlfile
      bamqc_zip: BamQC/zipfile
      readqc_zip: ReadQC/zipfile
      readqc_html: ReadQC/htmlfile
      macsDir: MACS-Auto/macsDir
      allmacsDir: MACS-All/macsDir
      nmmacsDir: MACS-NM/macsDir
      rpmwig: WIG-Auto/rpmwig
      outBW: WIG-Auto/outBW
      outtdf: WIG-Auto/outtdf
      allrpmwig: WIG-All/rpmwig
      alloutBW: WIG-All/outBW
      allouttdf: WIG-All/outtdf
      nmrpmwig: WIG-NM/rpmwig
      nmoutBW: WIG-NM/outBW
      nmouttdf: WIG-NM/outtdf
      bedfasta: MOTIFS/bedfasta
      flankbed: FlankBED/outfile
      memechipdir: MOTIFS/memechipdir
      summitmemechipdir: SummitMOTIFS/memechipdir
      amedir: MOTIFS/amedir
      summitamedir: SummitMOTIFS/amedir
      metagenesDir: MetaGene/metagenesDir
      sicerDir: SICER/sicerDir
      roseoutput: ROSE/RoseDir
      statsfile: PeaksQC/statsfile
      htmlfile: PeaksQC/htmlfile
      textfile: PeaksQC/textfile
    out: [ finalDir ]
    run: movealloutput.cwl
    scatter: [ sam_sort, fastq_metrics, rmdup_bam, bklist_bam, bamqc_html, bamqc_zip, readqc_zip, readqc_html, macsDir, allmacsDir, nmmacsDir, rpmwig, outBW, outtdf, allrpmwig, alloutBW, allouttdf, nmrpmwig, nmoutBW, nmouttdf, bedfasta, flankbed, memechipdir, summitmemechipdir, amedir, summitamedir, metagenesDir, sicerDir, roseoutput, statsfile, htmlfile, textfile ]
    scatterMethod: dotproduct
