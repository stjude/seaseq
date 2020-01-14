#!/usr/bin/perl
# process ChipSeq-1st-mapping output to second stage ChipSeq-2nd-peakcalls
# process outputfiles to specified folder

use strict;
use JSON;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
 
my ($help, $manual, $infile, $outfile, $json, $step, $folder);

GetOptions ("i|in=s"=>\$infile,"o|out=s"=>\$outfile, "s|step=s"=>\$step, "f|folder=s"=>\$folder);
my $usage = "perl $0 -i <log output file> -o <output file> -s <step> -f <folder>\n";
unless ( $infile && $step ) { die $usage; }

local $/; #Enable 'slurp' mode
open my $fh, "<", $infile;
$json = <$fh>;
close $fh;
my $data = decode_json($json);

if ($step == 1) {
  #output to tempfile: bamfile, zipfile, STATbkoutfile, STATbamoutfile, STATrmdupoutfile
  unless ($folder) { $folder = (split("\_fastq", $data->{'readqc_zip'}->{'nameroot'}))[0];  }
  open my $fh, ">>", $outfile;
  print $fh "\nbklistbamfile: \n  class: File\n  path: " .
      $data->{'bklist_bam'}->{'path'} . "\n";
  print $fh "\nbklistbamindex: \n  class: File\n  path: " .
      $data->{'bklist_index'}->{'path'} . "\n";
  print $fh "\nrmdupbamfile: \n  class: File\n  path: " .
      $data->{'rmdup_bam'}->{'path'} . "\n";
  print $fh "\nreadzipfile: \n  class: File\n  path: " .
      $data->{'readqc_zip'}->{'path'} . "\n";
  print $fh "\nSTATbkoutfile: \n  class: File\n  path: " .
      $data->{'stat_bk'}->{'path'} . "\n";
  print $fh "\nSTATbamoutfile: \n  class: File\n  path: " .
      $data->{'stat_bam'}->{'path'} . "\n";
  print $fh "\nSTATrmdupoutfile: \n  class: File\n  path: " .
      $data->{'stat_rmdup'}->{'path'} . "\n";
  print $fh "\nfastqmetricsfile: \n  class: File\n  path: " .
      $data->{'fastq_metrics'}->{'path'} . "\n";
  close $fh;

  my $newpath = (fileparse($data->{'bklist_bam'}->{'path'}))[1];
  #print "This is the path $newpath\n"; 
  `mkdir -p $folder/QC_files/STATS $folder/QC_files/FASTQC $folder/BAM_files`;

  #copy the relevant files to the specified folder
  `cp -rf $newpath/*fastqc* $folder/QC_files/FASTQC`;
  `cp -rf $newpath/*bam $newpath/*bai $folder/BAM_files`;
  `cp -rf $newpath/*metrics.txt $folder/QC_files/STATS`;
 
  print $folder;
}

elsif ($step == 2 && $folder) {
  #output to desired folder
  my $newpath = (fileparse($data->{'statsfile'}->{'path'}))[1];
  `mkdir -p $folder/PEAKS_files/NARROW_peaks $folder/PEAKS_files/BROAD_peaks $folder/PEAKS_files/STITCHED_REGIONS`;
  `mkdir -p $folder/PEAKSDisplay_files $folder/BAMDensity_files $folder/MOTIFS_files $folder/QC_files/STATS`;
  `cp -rf $newpath/ROSE_out/* $folder/PEAKS_files/STITCHED_REGIONS/`;
  `cp -rf $newpath/SICER_out/* $folder/PEAKS_files/BROAD_peaks`;
  `cp -rf $newpath/*ame_out $folder/MOTIFS_files`;
  `cp -rf $newpath/*memechip_out $folder/MOTIFS_files`;
  `cp -rf $newpath/bamdensity_out/* $folder/BAMDensity_files`;
  `cp -rf $newpath/*-stats* $folder/QC_files/STATS`;
  `cp -rf $newpath/*-p9_kd* $newpath/*-nm $folder/PEAKS_files/NARROW_peaks`;
  `cp -rf $newpath/*.wig* $newpath/*.bw $newpath/*.tdf $folder/PEAKSDisplay_files`;
}
