#!/usr/bin/bash

#getting the readLength of the fastq reads
fastqfile=$1;

readlength=$(zcat $fastqfile | head -n 2 | tail -n 1 | awk '{print length}');
echo $readlength  > $readlength.rdl;
