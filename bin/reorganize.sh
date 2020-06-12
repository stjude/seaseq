#!/bin/bash
#
# Reorganize output directory from ChipSEQ pipeline into preferred directory structure.
# ["final outdir"]
#

if [ $# -lt 2 ]; then
  echo ""
  echo 1>&2 Usage: $0 ["current scatter outdir"] ["final outdir"] 
  echo ""
  exit 1
fi


for folderpath in $1/*;
do
  newfolder=${folderpath##*/}
  mkdir -p $newfolder/QC_files/STATS $newfolder/QC_files/FASTQC $newfolder/BAM_files
  mkdir -p $newfolder/PEAKS_files/NARROW_peaks $newfolder/PEAKS_files/BROAD_peaks $newfolder/PEAKS_files/STITCHED_REGIONS
  mkdir -p $newfolder/PEAKSDisplay_files $newfolder/BAMDensity_files $newfolder/MOTIFS_files $newfolder/QC_files/STATS

  cp -rf $folderpath/*fastqc* $newfolder/QC_files/FASTQC
  cp -rf $folderpath/*bam $newfolder/BAM_files
  cp -rf $folderpath/*metrics.txt $newfolder/QC_files/STATS
  cp -rf $folderpath/ROSE_out/* $newfolder/PEAKS_files/STITCHED_REGIONS/
  cp -rf $folderpath/SICER_out/* $newfolder/PEAKS_files/BROAD_peaks
  cp -rf $folderpath/*ame_out $newfolder/MOTIFS_files
  cp -rf $folderpath/*memechip_out $newfolder/MOTIFS_files
  cp -rf $folderpath/bamdensity_out/* $newfolder/BAMDensity_files
  cp -rf $folderpath/*-stats* $newfolder/QC_files/STATS
  cp -rf $folderpath/*-p9_kd* $folderpath/*-nm $newfolder/PEAKS_files/NARROW_peaks
  cp -rf $folderpath/*.wig* $folderpath/*.bw $folderpath/*.tdf $newfolder/PEAKSDisplay_files

  mkdir -p $2
  mv $newfolder $2
done
