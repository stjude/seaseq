
# SEAseq and PEAseq "sub-sampled" dataset example run.
Please consult [documentation](https://github.com/stjude/seaseq/tree/master/docs#readme) for more details.

## Prerequisites
1. Docker (installed and running)
2. Cromwell.jar (> 51)

## Usage
### 1. Local OS (NOT recommended for real datasets)
``` sh
# replace <cromwell.jar> with actual cromwell.jar version downloaded
java -jar <cromwell.jar> run ../seaseq-case.wdl -i inputs-s_case.json -o options-s_case.json
```

### 2. HPC 

``` bash
# replace <hpc.config> with the applicable hpc configuration file
# replace <cromwell.jar> with actual cromwell.jar version downloaded
java -Dconfig.file=<hpc.config> -jar <cromwell.jar> run ../seaseq-case.wdl -i inputs-s_case.json -o options-s_case.json
```
View our provided list of [hpc configuration](https://github.com/stjude/seaseq/tree/master/docs/hpc-configurations) files.

## Further Details
Consult bash scripts (\*.sh) provided for more details of actual run.
