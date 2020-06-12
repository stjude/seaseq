#!/usr/bin/bash
# seaseq_wrap.sh
# ====== README
# This is the seaseq wrapper script written for
#   St. Jude, Abraham's Lab.
# This script works with St. Jude, HPC LSF cluster.
# ------
# To Run on St. Jude LSF.
#   bsub -P watcher -q compbio -J seaseq -o seaseq-out -e seaseq-err -N ./seaseq_wrap.sh <yml> <OUTPUTFOLDER>
#   bsub -P watcher -q compbio -J seaseq -o seaseq-out -e seaseq-err -N ./seaseq_wrap.sh inputparameters.yml OUTPUTFOLDER
# ------
# Programs and versions used
#   bowtie v. 1.2.2
#   fastqc v. 0.11.5
#   samtools v. 1.9
#   R v. 3.6.1
#   macs v. 041014
#   SICER2 v. 1.0.2
#   meme v. 4.11.2
#   spp v. 1.16.0
#   bedtools/2.25.0
#   python v. 3.7.0
#   java v. 1.8.0_60
#   perl v. 5.10.1
#   wigToBigWig v. 4
#   bedops v. 2.4.2
#   igvtools v. 2.3.2
#   ROSE v. 1.1.0
#   BAM2GFF v. 1.1.0
# ------
# Args:
#   seaseqroot : PATH to workflow location
#   inputyml : input yaml containing all required variables
#   outputfolder : final output directory
# ======

# ====== VARIABLES
# Required Args
# ------
if [ $# -lt 2 ]; then
  echo ""
  echo 1>&2 Usage: $0 ["YML"] ["OUTPUT FOLDER"]
  echo ""
  exit 1
fi
actualpath=$(realpath $0)
seaseqroot=${actualpath%/*}
inputyml=$1
outputfolder=$2
# ------
# Permanent Args
# ------
seaseqworkflow="$seaseqroot/cwl/seaseq_pipeline.cwl"
NEW_UUID=${NEW_UUID:=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 5 | head -n 1)"_"`date +%s`} #temporary file for the 2nd step

out="$(pwd)/seaseq-"$NEW_UUID"-outdir"
tmp="$(pwd)/seaseq-"$NEW_UUID"-tmpdir"
jobstore="seaseq-"$NEW_UUID"-jobstore"
logtxt="seaseq-"$NEW_UUID"-log.txt"
logout="seaseq-"$NEW_UUID"-log_out"
logerr="seaseq-"$NEW_UUID"-log_err"

export PATH=$PATH:$seaseqroot/bin
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

cwltool \
--no-container \
--preserve-entire-environment \
--copy-outputs \
--leave-tmpdir \
--tmpdir-prefix $tmp \
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
  OUTPUTFOLDER=$(chromatinSEreadjson.pl -i $logout -s 1)
  chromatinSEreadjson.pl -i $logout -s 2 -f $OUTPUTFOLDER

  #rm -rf *$NEW_UUID*
  mkdir -p $outputfolder
  mv $OUTPUTFOLDER $outputfolder

  echo "STATUS:  Cleaned Up All files with $NEW_UUID"
  echo "SUCCESS: CHIPSEQ - SE Pipeline Completed"
else
  echo "ERROR:   ChipSeq-ALL workflow terminated with errors"
fi

# ======
