#!/bin/bash
## Test WDL script using samplefiles provided.

# CROMWELL on local machine
local_cromwell="/Users/madetunj/Downloads/cromwell-61.jar"

# CROMWELL on St. Jude hpc user directory
sjhpc_cromwell="/home/madetunj/madetunj/software/cromwell-52.jar"
lsf_config="/home/madetunj/.commands/lsf.conf"
input="inputs-p_control.json"
option="options-p_control.json"

# STD OUT and ERR files
logout="wdlpeaseq-control_out"
logerr="wdlpeaseq-control_err"

# check if cromwell app exists
if [ -f "$local_cromwell" ]; then
    java -jar $local_cromwell \
        run ../peaseq-control.wdl \
        -i $input \
        -o $option \
        1>$logout 2>$logerr

elif [ -f "$sjhpc_cromwell" ]; then 
    #script syntax
    wdlscript="java -Dconfig.file=$lsf_config \
        -jar $sjhpc_cromwell \
        run ../peaseq-control.wdl \
        --inputs $input \
        --options $option"
    bsub -P watcher -q compbio \
        -R "rusage[mem=10000]" \
        -J wdlpeaseq-control \
        -o $logout \
        -e $logerr \
        -N $wdlscript

else
    echo "cromwell executable jar doesN'T exist"

fi

