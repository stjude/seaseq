#!/usr/bin/perl
# process ChipSeq-1st-mapping output to second stage ChipSeq-2nd-peakcalls
# process outputfiles to specified folder

use strict;
use JSON;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use Data::Dumper qw(Dumper);

my ($help, $manual, $infile, $outfile, $json, $step, $folder, $exec, $toil, $newpath);

GetOptions ("i|in=s"=>\$infile,"o|out=s"=>\$outfile, "s|step=s"=>\$step, "f|folder=s"=>\$folder, "exec"=>\$exec, "toil"=>\$toil);
my $usage = "perl $0 -i <log output file> -o <output file> -s <step> -f <folder> [-exec/-toil]\n";
unless ( $infile && $step ) { die $usage; }

#reading the json file
local $/; #Enable 'slurp' mode
open my $fh, "<", $infile;
$json = <$fh>;
close $fh;
my $data = decode_json($json);
#print Data::Dumper->Dump( [ \$data ], [ qw(*GENE) ] );

#After the 1st stage (1st-mapping) is completed
if ($step == 1) {
  #output to tempfile: bamfile, zipfile, STATbkoutfile, STATbamoutfile, STATrmdupoutfile
  unless ($folder) { $folder = (split("\_fastq", $data->{'readqc_zip'}->{'nameroot'}))[0];  }
  open my $fh, ">>", $outfile;

  if ($toil) {
    print $fh "\nbklistbamfile: \n  class: File\n  path: " .
        (split("file://", $data->{'bklist_bam'}->{'location'}))[1] . "\n";
    print $fh "\nrmdupbamfile: \n  class: File\n  path: " .
        (split("file://", $data->{'rmdup_bam'}->{'location'}))[1] . "\n";
    print $fh "\nreadzipfile: \n  class: File\n  path: " .
        (split("file://", $data->{'readqc_zip'}->{'location'}))[1] . "\n";
    print $fh "\nSTATbkoutfile: \n  class: File\n  path: " .
        (split("file://", $data->{'stat_bk'}->{'location'}))[1] . "\n";
    print $fh "\nSTATbamoutfile: \n  class: File\n  path: " .
        (split("file://", $data->{'stat_bam'}->{'location'}))[1] . "\n";
    print $fh "\nSTATrmdupoutfile: \n  class: File\n  path: " .
        (split("file://", $data->{'stat_rmdup'}->{'location'}))[1] . "\n";
    print $fh "\nfastqmetricsfile: \n  class: File\n  path: " .
        (split("file://", $data->{'fastq_metrics'}->{'location'}))[1] . "\n";
    close $fh;
    $newpath = (fileparse((split("file://", $data->{'bklist_bam'}->{'location'}))[1]))[1];
  }
  else { #if $exec or nothing is specified 
    print $fh "\nbklistbamfile: \n  class: File\n  path: " .
        $data->{'bklist_bam'}->{'path'} . "\n";
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
    $newpath = (fileparse($data->{'bklist_bam'}->{'path'}))[1];    
  }
  `mkdir -p $folder/QC_files/STATS $folder/QC_files/FASTQC $folder/BAM_files`;

  #copy the relevant files to the specified folder
  `cp -rf $newpath/*fastqc* $folder/QC_files/FASTQC`;
  `cp -rf $newpath/*bam $newpath/*bai $folder/BAM_files`;
  `cp -rf $newpath/*metrics.txt $folder/QC_files/STATS`;
 
  print $folder;
} #end if for step1

#After the 2nd stage (nd-peakcalls) is completed
elsif ($step == 2 && $folder) {
  #output to desired folder
  if ($toil) {
    $newpath = (fileparse((split("file://", $data->{'statsfile'}->{'location'}))[1]))[1];
  }
  else { #if $exec or nothing is specified 
    $newpath = (fileparse($data->{'statsfile'}->{'path'}))[1];
  }
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
} #end if for step2
