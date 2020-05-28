#!/usr/bin/bash
# seaseq_wrap.sh

# ====== README
# This is the seaseq wrapper with raptr and sqlite connections
#  script written for St. Jude, Abraham's Lab.
# This script works with St. Jude, HPC LSF cluster.
# This script require Toil and sjhpc modules.
# ------
# To Run on St. Jude LSF.
#   bsub -P watcher -q compbio -J seaseq -o seaseq-out -e seaseq-err -N ./seaseq_wrap.sh <yml> <OUTPUTFOLDER>
#   bsub -P watcher -q compbio -J seaseq -o seaseq-out -e seaseq-err -N ./seaseq_wrap.sh inputparameters.yml OUTPUTFOLDER
# ------
# Required variables:
#   SQLITE : PATH to sqlite database
#   seaseqroot : PATH to workflow location
#
# Command Line variables:
#   jobid : Job ID to update job status in sqlite
#   inputyml : input yaml containing all required variables
#   outputfolder : final output directory
# ======

# ====== VARIABLES
# Required Args
# ------
seaseqroot="/research/rgs01/project_space/abrahgrp/Software_Dev_Sandbox/common/madetunj/ChipSeqPipeline"
SQLITE="/research/rgs01/project_space/abrahgrp/FRAMEWORKOUTPUT/abrahgrp/.sqlite_seaseq.db"

# ======

if [ $# -lt 3 ]; then
  echo ""
  echo 1>&2 Usage: $0 ["JOB_ID"]["YML"] ["OUTPUT FOLDER"]
  echo ""
  exit 1
fi

# ====== UPDATE JOB STATUS IN SQLITE DATABASE
# Update job_status in sqlite database to 'processing'
# ------
module load sqlite/3.23.1

sqlite3 $SQLITE "UPDATE job_status SET js_status = 'processing' WHERE \
         job_id = $1;"
# ======
inputyml=$2
outputfolder=$3
# ------
# Permanent Args
# ------
seaseqworkflow="$seaseqroot/bin/seaseq_pipeline.cwl"
NEW_UUID=${NEW_UUID:=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 5 | head -n 1)"_"`date +%s`} #temporary file for the 2nd step

out="$(pwd)/seaseq-"$NEW_UUID"-outdir"
tmp="$(pwd)/seaseq-"$NEW_UUID"-tmpdir"
jobstore="seaseq-"$NEW_UUID"-jobstore"
logtxt="seaseq-"$NEW_UUID"-log.txt"
logout="seaseq-"$NEW_UUID"-log_out"
logerr="seaseq-"$NEW_UUID"-log_err"

export PATH=$PATH:$seaseqroot/scripts
# ------
# Loading St. Jude HPC cluster programs
# ------
source $seaseqroot/sjhpc_modules.ml

# ======

# ====== WORKFLOW EXECUTION
# Toil Workflow Execution on St. Jude HPC Cluster
# ------
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
$seaseqworkflow $inputyml 1>$logout 2>$logerr
# ------
# Extract and Clean Up output files.
#   Relevant files are moved to indicated "$outputfolder"
#   Remove all temporary files.
# ------
if [ -s $logout ]
then

  ##extracting relevant files from 1st step to the next step & to outputfolder
  OUTPUTFOLDER=$(chromatinSEreadjson.pl -i $logout -s 1 -toil)
  chromatinSEreadjson.pl -i $logout -s 2 -f $OUTPUTFOLDER -toil

  #clean temporary files & move files to specified folder
  rm -rf *$NEW_UUID*
  mkdir -p $outputfolder
  mv $OUTPUTFOLDER $outputfolder

  #update job status and sqlite job_status
  echo "STATUS:  Cleaned Up All files with $NEW_UUID"
  echo "SUCCESS: CHIPSEQ - SE Pipeline Completed"
  sqlite3 $SQLITE "UPDATE job_status SET js_status = 'done' WHERE \
           job_id = $1;"
else
  echo "ERROR:   ChipSeq-ALL workflow terminated with errors"
  sqlite3 $SQLITE "UPDATE job_status SET js_status = 'failed' WHERE \
           job_id = $1;"
fi

# ======
