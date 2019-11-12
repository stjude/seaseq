#!/usr/bin/perl
# process ChipSeq-1st-mapping output to second stage ChipSeq-2nd-peakcalls
# process outputfiles to specified folder

use strict;
use JSON;
use Getopt::Long;
use Pod::Usage;
 
my ($help, $manual, $infile, $outfile, $json, $step, $folder);

GetOptions ("i|in=s"=>\$infile,"o|out=s"=>\$outfile, "s|step=s"=>\$step, "f|folder=s"=>\$folder);
my $usage = "perl $0 -i <log output file> -o <output file> -s <step>\n";
unless ( $infile && $step ) { die $usage; }

local $/; #Enable 'slurp' mode
open my $fh, "<", $infile;
$json = <$fh>;
close $fh;
my $data = decode_json($json);

if ($step == 1) {
  #output to tempfile: bamfile, zipfile, STATbkoutfile, STATbamoutfile, STATrmdupoutfile
  open my $fh, ">>", $outfile;
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
  close $fh;

  if ($folder) {
    `mkdir -p $folder`;
    #copy the relevant files to the specified folder
    `cp $data->{'sam_sort'}->{'path'} $folder`;
    `cp $data->{'rmdup_bam'}->{'path'} $folder`;
    `cp $data->{'rmdup_index'}->{'path'} $folder`;
    `cp $data->{'bklist_bam'}->{'path'} $folder`;
    `cp $data->{'bklist_index'}->{'path'} $folder`;
    `cp $data->{'readqc_zip'}->{'path'} $folder`;
    `cp $data->{'bamqc_zip'}->{'path'} $folder`;
    `cp $data->{'readqc_html'}->{'path'} $folder`;
    `cp $data->{'bamqc_html'}->{'path'} $folder`;
  }
}

elsif ($step == 2 && $folder) {
  #output to desired folder
  `cp $data->{'sicerbed'}->{'path'} $folder`;
  `cp $data->{'sicergraph'}->{'path'} $folder`;
  `cp $data->{'sicerisland'}->{'path'} $folder`;
  `cp $data->{'genebodypdf'}->{'path'} $folder`;
  `cp $data->{'promoterspdf'}->{'path'} $folder`;
  `cp $data->{'genebodyheatmappng'}->{'path'} $folder`;
  `cp $data->{'promotersheatmappng'}->{'path'} $folder`;
  
  #auto
  `cp $data->{'summits'}->{'path'} $folder`;
  `cp $data->{'ameout'}->{'path'} $folder`;
  `cp $data->{'memeout'}->{'path'} $folder`;
  `cp $data->{'memehtml'}->{'path'} $folder`;
  `cp $data->{'amehtml'}->{'path'} $folder`;
  `cp $data->{'outBW'}->{'path'} $folder`;
  `cp $data->{'peaksxls'}->{'path'} $folder`;
  `cp $data->{'peaksbed'}->{'path'} $folder`;
  `cp $data->{'treatwig'}->{'path'} $folder`;
  `cp $data->{'statout'}->{'path'} $folder`;
  `cp $data->{'outtdf'}->{'path'} $folder`;
  `cp $data->{'summitameout'}->{'path'} $folder`;
  `cp $data->{'summitmemeout'}->{'path'} $folder`;
  `cp $data->{'summitmemehtml'}->{'path'} $folder`;
  `cp $data->{'summitamehtml'}->{'path'} $folder`;

  #all
  `cp $data->{'allsummits'}->{'path'} $folder`;
  `cp $data->{'allameout'}->{'path'} $folder`;
  `cp $data->{'allmemeout'}->{'path'} $folder`;
  `cp $data->{'allmemehtml'}->{'path'} $folder`;
  `cp $data->{'allamehtml'}->{'path'} $folder`;
  `cp $data->{'alloutBW'}->{'path'} $folder`;
  `cp $data->{'allpeaksxls'}->{'path'} $folder`;
  `cp $data->{'allpeaksbed'}->{'path'} $folder`;
  `cp $data->{'alltreatwig'}->{'path'} $folder`;
  `cp $data->{'allstatout'}->{'path'} $folder`;
  `cp $data->{'allouttdf'}->{'path'} $folder`;
  `cp $data->{'allsummitameout'}->{'path'} $folder`;
  `cp $data->{'allsummitmemeout'}->{'path'} $folder`;
  `cp $data->{'allsummitmemehtml'}->{'path'} $folder`;
  `cp $data->{'allsummitamehtml'}->{'path'} $folder`;
}
