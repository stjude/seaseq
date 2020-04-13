# Chromatin SE analysis pipeline.
- Currently designed as ChipSeq

# (S)ingle (E)nd (A)ntibody (SEQ)uencing pipeline

## USAGE:

Three workflows are available to run SEASEQ pipeline.

1. BAM mapping & Peak Calling.

1. All in one pipeline.

1. All in one pipeline for multiple fastQ.


## To RUN:

1. The two-step pipeline (BAM mapping & Peak Calling):

	1. (CWLExec) JobSubmit.sh
	
		``` 
		bsub -R -P watcher -q compbio -J exec-cwl -o exec-cwl_out -e exec-cwl_err -N ./JobSubmit.sh
		```

	1. (Toil) ToilJob.sh
		
		```
		bsub -P watcher -q compbio -J toil-cwl -o toil-cwl_out -e toil-cwl_err -N ./ToilJob.sh
		```
1. The single-step pipeline:
	
	1. (CWLExec) ALLExecJob.sh
	
		```
		bsub -P watcher -q compbio -J allexec -o allexec_out -e allexec_err -N ./ALLExecJob.sh
		```

	1. (Toil) ALLToilJob.sh
	
		```
		bsub -P watcher -q compbio -J alltoil -o alltoil_out -e alltoil_err -N ./ALLToilJob.sh
		```
		
1. The single-step pipeline for multiple FastQs.
	
	1. (Toil) ScatterToil.sh
	
		```
		bsub -P watcher -q compbio -J scattertoil -o scattertoil_out -e scattertoil_err -N ./ScatterToil.sh [input YAML File] [OUTPUTFOLDER]
		```

