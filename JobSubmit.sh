#!/usr/bin/env bash
#------
###SYNTAX to run
#bsub -R "rusage[mem=10000]" -P watcher -q compbio -J exec-cwl -o exec-cwl_out -e exec-cwl_err -N ./JobSubmit.sh
####

#------
###FILES
#------
#location="/rgs01/project_space/abrahgrp/Software_Dev_Sandbox/common/madetunj/ChipSeqPipeline"
location="/rgs01/project_space/abrahgrp/Software_Dev_Sandbox/common/madetunj/SEQ2"
parameters="$location/inputparameters.yml"
config="$location/LSFconfig.json"
firstscript="$location/workflows/ChromatinSE-1st-mapping.cwl"
secondscript="$location/workflows/ChromatinSE-2nd-peakcalls.cwl"

OLD_UUID=$1
NEW_UUID=${NEW_UUID:=${OLD_UUID%%.yml}} #reuse old file
NEW_UUID=${NEW_UUID:=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 5 | head -n 1)"_"`date +%s`} #temporary file for the 2nd step

#temporary output & error files
out="$(pwd)/exec-"$NEW_UUID"-outdir"
tmp="$(pwd)/exec-"$NEW_UUID"-tmpdir"
logout="chromatinSE-"$NEW_UUID"-log_out"
logerr="chromatinSE-"$NEW_UUID"-log_err"

#------
###Modules & PATH update
#------
module load /rgs01/project_space/zhanggrp/MethodDevelopment/common/CWL/modulefiles/cwlexec/latest
module load node
module load igvtools/2.3.2
module load fastqc/0.11.5
module load bowtie/1.2.2
module load samtools/1.9
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

export PATH=$PATH:$location/scripts
export R_LIBS_USER=$R_LIBS_USER:/rgs01/project_space/abrahgrp/Software_Dev_Sandbox/common/madetunj/R #to find SPP local package

#------
###WORKFLOW
#------
##cwlexec 1st step
echo "STATUS:  Temporary files named with $NEW_UUID"
echo "UPDATE:  STEP1 in progress"

if [ ! -f $NEW_UUID.yml ] 
then
  ##removing work and out files
  mkdir -p $tmp $out
  
  ##excuting step one
  cwlexec -p -w $tmp -o $out -c $config -p $firstscript $parameters 1>$logout 2>$logerr
else
  echo "NOTICE:  STEP1 initially completed with temp yml $NEW_UUID.yml"
fi

##extract 1st section
if [ -s $logout ]
then
  echo "UPDATE:  STEP1 Successfully completed"
  cp -f $parameters $NEW_UUID.yml

  ##extracting relevant files from 1st step to the next step & to outputfolder
  OUTPUTFOLDER=$(chromatinSEreadjson.pl -i $logout -o $NEW_UUID.yml -s 1)

  echo "UPDATE:  STEP2 in progress"


  ##cwlexec 2nd step
  cwlexec -p -w $tmp -o $out -c $config -p $secondscript $NEW_UUID.yml 1>$logout.2 2>$logerr.2

  ##extract require files into folder
  if [ -s $logout.2 ]
  then
    echo "UPDATE:  STEP2 Successfully completed"
    chromatinSEreadjson.pl -i $logout.2 -s 2 -f $OUTPUTFOLDER

    #log and error files are saved
    mkdir -p $OUTPUTFOLDER/Log_files
    mv $NEW_UUID.yml $logerr $logerr.2 $logout $logout.2 $OUTPUTFOLDER/Log_files

    echo "UPDATE:  CHIPSEQ - SE Pipeline Completed"

  else
    echo "ERROR:   STEP2 for ChipSeq workflow terminated with errors"
  fi
else
  echo "ERROR:   STEP1 for ChipSeq workflow terminated with errors"
fi

