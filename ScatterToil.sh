#!/usr/bin/env bash
#------
###SYNTAX to run
# bsub -P watcher -q compbio -J scattertoil -o scattertoil_out -e scattertoil_err -N ./ScatterToil.sh <yml> <OUTPUTFOLDER>
# example:
#  bsub -P watcher -q compbio -J scattertoil -o scattertoil_out -e scattertoil_err -N ./ScatterToil.sh inputparameters-scatter.yml OUTPUTFOLDER

location="/rgs01/project_space/abrahgrp/Software_Dev_Sandbox/common/madetunj/ChipSeqPipeline"
script="$location/workflows/ChromatinSE-Scatter.cwl"

NEW_UUID=${NEW_UUID:=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 5 | head -n 1)"_"`date +%s`} #temporary file for the 2nd step

#------
###FILES
#------
if [ $# -lt 2 ]; then
  echo ""
  echo 1>&2 Usage: $0 ["YML"] ["OUTPUT FOLDER"] 
  echo ""
  exit 1
fi

parameters=$1
OUTPUTFOLDER=$2

#temporary output & error files
out="$(pwd)/scatter-"$NEW_UUID"-outdir"
tmp="$(pwd)/scatter-"$NEW_UUID"-tmpdir"
jobstore="scatter-"$NEW_UUID"-jobstore"
logtxt="scatter-"$NEW_UUID"-log.txt"
logout="scatter-"$NEW_UUID"-log_out"
logerr="scatter-"$NEW_UUID"-log_err"

#------
###Modules & PATH update
#------
module load node
module load igvtools/2.3.2
module load fastqc/0.11.5
module load bowtie/1.2.2
module load macs/041014
module load ucsc/041619
module load R/3.6.1
module load bedtools/2.25.0
module load meme/4.11.2
module load bedops/2.4.2
module load java/1.8.0_60
module load BAM2GFF/1.1.0
module load ROSE/1.1.0
module load SICER2/1.0.1
module load samtools/1.9

export PATH=$PATH:$location/scripts
export R_LIBS_USER=$R_LIBS_USER:/rgs01/project_space/abrahgrp/Software_Dev_Sandbox/common/madetunj/R #local SPP

#------
###WORKFLOW
#------
echo "STATUS:  Temporary files named with $NEW_UUID"
mkdir -p $tmp $out
 
rm -rf $jobstore $logtxt

toil-cwl-runner --batchSystem=lsf \
--preserve-entire-environment \
--disableCaching \
--logFile $logtxt \
--jobStore $jobstore \
--clean never \
--workDir $tmp \
--cleanWorkDir never \
--outdir $out \
$script $parameters 1>$logout 2>$logerr

##extract 1st section
if [ -s $logout ]
then

  # moving relevant file to desired output directory
  reorganize.sh $out $OUTPUTFOLDER
  rm -rf *$NEW_UUID*

  echo "Results stored in $OUTPUTFOLDER"
  echo "SUCCESS:  CHIPSEQ - SE Pipeline Completed"

else
  echo "ERROR:   ChipSeq-ALL workflow terminated with errors"
fi

