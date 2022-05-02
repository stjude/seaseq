# HPC Cluster Schedulers

HPC Systems using Singularity will require a cromwell configuration file. 
The configuration file contains the needed commands to submit and monitor singularity executions. 

We provided configuration files for the different schedulers.

1. IBM SPECTRUM LSF
	i. `lsf.conf` : The lsf.conf file should be included for every run.
	```bash
	# For sample FASTQs only
	java -Dconfig.file=lsf.conf -jar cromwell.jar run seaseq-case.wdl -i inputs.json -o options.json
	
	# For sample FASTQs with Input Control or IgG
	java  -Dconfig.file=lsf.conf -jar cromwell.jar run seaseq-control.wdl -i inputs.json -o options.json
	```

	ii. `check-job-alive` : This script should be included in your `$PATH`
	```bash
	# Run the following commands to put it in your $PATH:
	wget https://raw.githubusercontent.com/stjude/seaseq/master/docs/hpc-configurations/check-job-alive
	mv check-job-alive ~/bin/
	```

1. SLURM
	i. `slurm-conf` : The slurm.conf file should be included for every run
        ```bash
        # For sample FASTQs only
        java -Dconfig.file=slurm.conf -jar cromwell.jar run seaseq-case.wdl -i inputs.json -o options.json
        
        # For sample FASTQs with Input Control or IgG
        java  -Dconfig.file=slurm.conf -jar cromwell.jar run seaseq-control.wdl -i inputs.json -o options.json
        ```
       


## Questions

Please contact Maintainer for further assistance
