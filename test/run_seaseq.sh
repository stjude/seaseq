#!/bin/bash
## Test WDL script using samplefiles provided.

# CROMWELL on local machine
local_cromwell="/Users/madetunj/Downloads/cromwell-51.jar"

# CROMWELL on St. Jude hpc user directory
sjhpc_cromwell="/home/madetunj/.software/cromwell-52.jar"
lsf_config="/home/madetunj/.commands/lsf.conf"

# STD OUT and ERR files
logout="wdlseaseq-out"
logerr="wdlseaseq-err"

# removing old files
rm -rf SEASEQ seaseq_logs wdlseaseq*

# check if cromwell app exists
if [ -f "$local_cromwell" ]; then
    java -jar $local_cromwell \
        run ../seaseq.wdl \
        -i seaseqinputs.json \
        -o seaseqoptions.json \
        1>$logout 2>$logerr

elif [ -f "$sjhpc_cromwell" ]; then 
    #script syntax
    wdlscript="java -Dconfig.file=$lsf_config \
        -jar $sjhpc_cromwell \
        run ../seaseq.wdl \
        --inputs seaseqinputs.json \
        --options seaseqoptions.json"
    bsub -P watcher -q compbio \
        -J wdlseaseq \
        -o $logout \
        -e $logerr \
        -N $wdlscript

else
    echo "cromwell executable jar doesN'T exist"

fi

