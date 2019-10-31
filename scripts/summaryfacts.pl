#!/usr/bin/perl
#summary stats of all files provided

use Pod::Usage;
use strict; 
use warnings;
use File::Basename;
use Getopt::Long;

my $PATH = "/rgs01/project_space/abrahgrp/Software_Dev_Sandbox/common/madetunj/software";
my ($help, $manual, $rmdupbam, $peaksbed, $bamflag, $rmdupflag, $bkflag, $outfile, $fastqczip);
my $usage = "perl $0 -b <bam file> -p <peaks bed> -bamflag <bamflagstat> -rmdupflag <rmdupflagstat> -bkflag <bklistflagstat> -outfile <outputfile> -fqc <fastqczipfile>\n";

GetOptions ("b|bam=s"=>\$rmdupbam,"p|peak=s"=>\$peaksbed,"bamflag=s"=>\$bamflag,"rmdupflag=s"=>\$rmdupflag,"bkflag=s"=>\$bkflag,"fqc|fastqczip=s"=>\$fastqczip,"outfile|o=s"=>\$outfile);
unless ($rmdupbam || $peaksbed || $bamflag || $rmdupflag || $bkflag || $fastqczip) { die $usage;}

#Filenames
my ($countsfile, $bambed, $sppout, $statsout);
my ($totalreads, $mappedreads, $rmdupreads, $bklistreads);

#Initialize variables
my ($Uniquecnt, $Totalcnt, $Fripcnt, $FRIP, $peaks, $NRF) = (0,0,0,0,0,0);
my $prev = "NA";

if ($outfile) { $statsout = fileparse($outfile, qr/\.[^.]*(\..*)?$/)."-stats.out"; } 
else { $statsout = fileparse($rmdupbam, qr/\.[^.]*(\..*)?$/)."-stats.out"; }

open (OUT, ">$statsout"); #opening outputfile

#working with FastQC
if ($fastqczip) {
  $totalreads = `unzip -p $fastqczip */fastqc_data.txt | grep "Total Sequences" | awk -F' ' '{print \$3}'`;
  print OUT "Total Reads = $totalreads";
}

#working with Flagstat
if ($bamflag) {
  $mappedreads = `head -n 5 $bamflag | tail -n 1 | awk -F" " '{print \$1}'`;
  print OUT "Total Mapped Reads = $mappedreads";
}

if ($bkflag) {
  $bklistreads = `head -n 5 $bkflag | tail -n 1 | awk -F" " '{print \$1}'`;
  print OUT "Total After RM BlackList regions = ",$bklistreads;
}

if ($rmdupflag) {
  $rmdupreads = `head -n 5 $rmdupflag | tail -n 1 | awk -F" " '{print \$1}'`;
  print OUT "Total After RM Duplicates = ",$rmdupreads;
}

#PROCESS NRF
if ($rmdupbam) {
  print "Processing NRF score ...";
  unless ($rmdupbam =~ /\.bam$/) { die $!; }
  $bambed = fileparse($rmdupbam, qr/\.[^.]*(\.bam)?$/)."-bam2bed.bed";
  $sppout = fileparse($rmdupbam, qr/\.[^.]*(\.bam)?$/)."-spp.out";

  `bamToBed -i $rmdupbam > $bambed`;
  open(IN, "<$bambed") || die $!;
  while(<IN>){
    chomp;
    my @line = split("\t",$_); 
    $Totalcnt++;
    my $t = join("_",@line[0..2]);
    $Uniquecnt++ unless($t eq $prev);
    $prev=$t;
  } close (IN);
  $NRF = $Uniquecnt/$Totalcnt;
  print OUT "From SPP: Mapped Reads = $Totalcnt\nFrom SPP: Unique Reads = $Uniquecnt\nNRF score = $NRF\n";
  print ".. Done\n";


#Process NSC (normalized) + RSC (relative strand cross-correlation coefficient)
  print "Processing NSC and RSC ...";
  print "run_spp.R -c=$rmdupbam -savp -out=$sppout 1>outfile.fake 2>outfile.fake;\n";
  `Rscript $PATH/run_spp.R -c=$rmdupbam -savp -out=$sppout 1>outfile.fake 2>outfile.fake;`;
  open (IN, "<$sppout"); # || die $!; 
  my ($NSC, $RSC) = (split("\t", <IN>))[8,9]; close(IN);
  print OUT "Normalized Strand cross-correlation coefficient (NSC) = $NSC\nRelative Strand cross-correlation Coefficient (RSC) = $RSC\n";
  print ".. Done\n";
}

if ($peaksbed && $rmdupbam) {
#FRIP score
print "Processing FRIP score ... ";
  $countsfile = fileparse($peaksbed, qr/\.[^.]*(\.bed)?$/)."-out.txt"; 
  `sort-bed $bambed | intersectBed -sorted -a $peaksbed -b stdin -c > $countsfile`;

  open(IN,"<$countsfile")|| die $!;
  while(<IN>){
    chomp;
    my @l=split("\t",$_);
    $Fripcnt+=$l[-1];
    $peaks++; 
  } close (IN); `rm -rf $countsfile`;
  $FRIP = sprintf ("%.6f", ($Fripcnt/$Totalcnt));
  print OUT "Total Peaks = $peaks\nFRIP score = $FRIP\n";

  print ".. Done\n";
  close OUT;
}

`rm -rf outfile.fake $sppout`;
