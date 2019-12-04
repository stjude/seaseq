#!/usr/bin/bash

#extracting mean from basic metrics file

datafile=$1;

taglength=$(tail -n 1 $datafile | awk '{print $4}');
echo $taglength  > $taglength.rdl;
