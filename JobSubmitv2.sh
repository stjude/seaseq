#!/usr/bin/env bash
#------
###SYNTAX to run
#bsub -P watcher -q compbio -R "rusage[mem=10000]" -J execv2-cwl -o execv2-cwl_out -e execv2-cwl_err -N ./JobSubmitv2.sh
####

#------
###FILES
#------
location="/rgs01/project_space/abrahgrp/Software_Dev_Sandbox/common/madetunj/ChipSeqPipeline"
config="$location/LSFconfig.json"
parameters="$location/inputparameters.yml"
firstscript="$location/workflows/ChipSeq-1st-mapping.cwl"
secondscript="$location/workflows/ChipSeq-2nd-peakcalls.cwl"
OUTPUTFOLDER="ChipSeq-SE_"`date +%Y%m%d-%H%M` #final outputfolder
NEW_UUID=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 13 | head -n 1) #temporary file for the 2nd step

#temporary output & error files
out="$(pwd)/outdir"
tmp="$(pwd)/tmpdir"
logout="chipseqoutfile_out"
logerr="chipseqoutfile_err"

#------
###Modules & PATH update
#------
module load /rgs01/project_space/zhanggrp/MethodDevelopment/common/CWL/modulefiles/cwlexec/latest
module unload python #unloading all python versions
module load igvtools/2.3.2
module load fastqc/0.11.5
module load python/2.7.2
module load bowtie/1.2.2
module load samtools/1.9
module load macs/041014
module load wigToBigWig/4
module load R/3.6.1
module load bedtools/2.25.0
module load meme/4.11.2
module load sambamba/0.5.5
module load bedops/2.4.2
module load java/1.8.0_60
 
export PATH=$PATH:$location/scripts

#------
###HOUSEKEEPING
#------
#removing work and out files
rm -rf $tmp $out
mkdir -p $tmp $out

#------
###WORKFLOW
#------
##cwlexec 1st step
echo "UPDATE:  STEP1 in progress"

cwlexec -p -w $tmp -o $out -c $config -p $firstscript $parameters 1>$logout 2>$logerr

##extract 1st section
if [ -s $logout ]
then
  cp -f $parameters $NEW_UUID.yml

  ##extracting relevant files from 1st step to the next step & to outputfolder
  chipseqreadjson.pl -i $logout -o $NEW_UUID.yml -s 1 -f $OUTPUTFOLDER

  echo "UPDATE:  STEP2 in progress";

  ##cwlexec 2nd step
  cwlexec -p -w $tmp -o $out -c $config -p $secondscript $NEW_UUID.yml 1>$logout.2 2>$logerr.2

  ##extract require files into folder
  if [ -s $logout.2]
  then
    chipseqreadjson.pl -i $logout.2 -s 2 -f $OUTPUTFOLDER
    rm -rf $NEW_UUID.yml $logerr $logerr.2 $logout $logout.2
  else
    echo "ERROR:   STEP2 for ChipSeq workflow terminated with errors"
  fi
else
  echo "ERROR:   STEP1 for ChipSeq workflow terminated with errors";
fi

