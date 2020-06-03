#!/usr/bin/bash
#
#getting basic metrics of fastq reads. 
#

if [ $# -lt 1 ]; then
  echo ""
  echo 1>&2 Usage: $0 ["FASTQ file"] 
  echo ""
  exit 1
fi

#File handling
fastqfile=$1
outputfile=$2

# Processing fastqfile
zcat $fastqfile | awk 'NR%4==2' | awk '{print length}' | sort -n > values.dat

#standard deviation
stddev=$(awk '{x+=$0;y+=$0^2}END{print sqrt(y/NR-(x/NR)^2)}' values.dat)
#awk '{x+=$0;y+=$0^2}END{print sqrt(y/NR-(x/NR)^2)}' values.dat) #population std dev. 

median=$(awk '{ a[i++]=$1; } END { print a[int(i/2)]; }' values.dat)

average=$(awk '{ sum += $1 } END { if (NR > 0) print sum / NR }' values.dat)

minimum=$(head -n 1 values.dat)

maximum=$(tail -n 1 values.dat)

Q1=$(awk 'BEGIN{c=0} {total[c]=$1; c++;} END{print total[int(NR*0.25 - 0.5)]}' values.dat)

Q3=$(awk 'BEGIN{c=0} {total[c]=$1; c++;} END{print total[int(NR*0.75 - 0.5)]}' values.dat)

IQR=$(echo "$Q3-$Q1" | bc)


echo Min.$'\t'1st Qu.$'\t'Median$'\t'Mean$'\t'3rd Qu.$'\t'Max.$'\t'StdDev.$'\t'IQR > $outputfile
echo $minimum$'\t'$Q1$'\t'$median$'\t'$average$'\t'$Q3$'\t'$maximum$'\t'$stddev$'\t'$IQR >> $outputfile

rm -rf values.dat

#echo "std deviation " $stddev
#echo "median " $median
#echo "average " $average
#echo "minimum " $minimum
#echo "maximum " $maximum
#echo "Q1 " $Q1
#echo "Q3 " $Q3
#echo "IQR " $IQR

