#!/usr/bin/env bash
#------
###SYNTAX to run
#bsub -P watcher -q compbio -R "rusage[mem=10000]" -J exec-cwl -o exec-cwl_out -e exec-cwl_err -N ./JobSubmit.sh
####

#------
###FILES
#------
location="/rgs01/project_space/abrahgrp/Software_Dev_Sandbox/common/ChipSeqPipeline"
config="$location/LSFconfig.json"
parameters="$location/inputparameters.yml"
script="$location/workflows/ChipSeq-workflow.cwl"


#------
###Modules & PATH update
#------
module load /rgs01/project_space/zhanggrp/MethodDevelopment/common/CWL/modulefiles/cwlexec/latest
module load fastqc/0.11.5
module load python/2.7.2
module load bowtie/1.2.2
module load samtools/1.9
module load macs/041014
module load wigToBigWig/4
module load bedtools/2.25.0
module load meme/4.11.2
module load java/1.8.0_60
 
export PATH=$PATH:$location/scripts

#------
###HOUSEKEEPING
#------
#removing work and out files
rm -rf $(pwd)/tmpdir $(pwd)/outdir
mkdir $(pwd)/tmpdir $(pwd)/outdir

#------
###WORKFLOW
#------

##cwlexec
cwlexec -p -w $(pwd)/tmpdir -o $(pwd)/outdir \
-c $config -p $script $parameters 1>zoutfile_out 2>zoutfile_err
