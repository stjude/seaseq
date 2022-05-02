# HPC Cluster Schedulers

HPC Systems using Singularity will require a cromwell configuration file to be included for every execution.

```bash
#config='lsf.conf' or 'slurm.conf'
java -Dconfig.file=$config -jar cromwell.jar run seaseq-case.wdl -i inputs.json -o options.json
```

The configuration file contains the needed commands to submit and monitor singularity executions. 

We provided configuration files for the different schedulers:

1. **IBM SPECTRUM LSF** requires the `lsf.conf` file, and a `check-job-alive` script to be included in your `$PATH`

	```bash
	# Run the following commands to put it in your $PATH:
	wget https://raw.githubusercontent.com/stjude/seaseq/master/docs/hpc-configurations/check-job-alive
	mv check-job-alive ~/bin/
	```

1. **SLURM** only requires the `slurm.conf` file.
