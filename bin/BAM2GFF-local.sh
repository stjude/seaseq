#!/bin/bash
#
# BAM TO GFF to detect Read density across gene transcripting regions
#
# Edited from Version 1 12/11/2019

##############################################################
# ##### Please replace PATHTO with your own directory ###### #
##############################################################
#PATHTO=/rgs01/project_space/abrahgrp/Software_Dev_Sandbox/common/madetunj/software/NEW/BAM2GFF
#PYTHONPATH=$PATHTO/lib
#export PYTHONPATH
#export PATH=$PATH:$PATHTO/bin

if [ $# -lt 4 ]; then
  echo ""
  echo 1>&2 Usage: $0 ["GTF file"] ["feature type"] ["BAM file"] ["CHROM SIZES"] ["SAMPLENAME"]
  echo ""
  exit 1
fi

#================================================================================
#Parameters for running

# GTF files
GTFFILE=$1

#FEATURE TYPE
FEATURE=$2
FEATURE=${FEATURE:=gene}

# BAM file
BAMFILE=$3

# CHROM SIZES
CHROMSIZES=$4

#sample name
SAMPLENAME=$5
TEMPFILE=${BAMFILE##*/}
SAMPLENAME=${SAMPLENAME:=${TEMPFILE%%.*}}

echo "#############################################"
echo "######            BAM2GFF v1           ######"
echo "#############################################"

echo "BAM file: $BAMFILE"
echo "FEATURE type: $FEATURE"
echo "Sample Name: $SAMPLENAME"
#================================================================================
#
# GENERATING GFF files for each genomic region
#
mkdir -p annotation
echo "Extracting GENE regions"
echo "BAM2GFF_gtftogenes.py -g $GTFFILE -f $FEATURE -c $CHROMSIZES"
BAM2GFF_gtftogenes.py -g $GTFFILE -f $FEATURE -c $CHROMSIZES
echo
#
# BAM TO GFF main code
#
mkdir -p matrix
echo "Working on Promoter Region"
echo "BAM2GFF_main.py -b $BAMFILE -i annotation/promoters.gff -m 100 -o matrix/promoters.txt"
BAM2GFF_main.py -b $BAMFILE -i annotation/promoters.gff -m 100 -o matrix/promoters.txt
echo
echo "Working on GeneBody Region"
echo "BAM2GFF_main.py -b $BAMFILE -i annotation/genes.gff -m 100 -o matrix/genebody.txt"
BAM2GFF_main.py -b $BAMFILE -i annotation/genes.gff -m 100 -o matrix/genebody.txt
echo
echo "Working on Upstream Region"
echo "BAM2GFF_main.py -b $BAMFILE -i annotation/upstream.gff -m 50 -o matrix/upstream.txt"
BAM2GFF_main.py -b $BAMFILE -i annotation/upstream.gff -m 50 -o matrix/upstream.txt
echo
echo "Working on Downstream Region"
echo "BAM2GFF_main.py -b $BAMFILE -i annotation/downstream.gff -m 50 -o matrix/downstream.txt"
BAM2GFF_main.py -b $BAMFILE -i annotation/downstream.gff -m 50 -o matrix/downstream.txt
echo
#
# PLOTS
#
echo "BAM2GFF_plots.R $SAMPLENAME"
BAM2GFF_plots.R $SAMPLENAME

echo "Done!"
