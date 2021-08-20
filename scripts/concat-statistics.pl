#!/usr/bin/perl
#summary stats of all files provided

use Pod::Usage;
use strict; 
use warnings;
use File::Basename;
use Getopt::Long;

my ($help, $manual, $sampleqc, $controlqc, $overallqc, $outfile);

my $usage = "perl $0 -s <sample config ml> -c <control config ml> -q <final seaseq config ml> -outfile <outputfile>\n";

GetOptions (
            "s|sample=s"=>\$sampleqc, "c|control=s"=>\$controlqc, "q|overall=s"=>\$overallqc,
            "outfile|o=s"=>\$outfile);

#Filenames
my ($htmlfile, $textfile);


#Initialize variables
my (%SQC, %CQC, %OQC);
my ($statsout);

#QC dictionary
my %stats = (
  1 => 'Overall_Quality',
  2 => 'Raw_Reads',
  3 => 'Base_Quality',
  4 => 'Sequence_Diversity',
  5 => 'Aligned_Percent',
  6 => 'Estimated_Fragment_Width',
  7 => 'Estimated_Tag_Length',
  8 => 'NRF',
  9 => 'PBC',
  10 => 'NSC',
  11 => 'RSC',
  12 => 'FRiP',
  13 => 'Total_Peaks',
  14 => 'Linear_Stitched_Peaks',
  15 => 'SE-like_Enriched_Regions');

#color names
my %color = ( "-2" => "#FF0000", "-1" => "#FF8C00", "0" => "#FFFF00", "1" => "#ADFF2F", "2" => "#008000" ); #red #orangered #yellow #greenyellow #green

#output file name
unless ($outfile) { 
  $statsout = "summarystats-stats.html"; 
} else {
  unless ($outfile =~ /\-stats.html$/) { $statsout = fileparse($outfile, qr/(\.[\w\d]+)$/)."-stats.html"; }
  else { $statsout = $outfile; }
}
#html output file name
$htmlfile = $statsout; $htmlfile =~ s/stats.html/stats.htmlx/; #creating html file
$textfile = $statsout; $textfile =~ s/stats.html/stats.txt/; #creating html file


#Reading Sample QC
open(SAMPLE, "<$sampleqc") || die $!;
while(<SAMPLE>){
  chomp $_;
  my @data = split("\t", $_);
  $SQC{$data[0]}{'value'} = $data[1];
  $SQC{$data[0]}{'score'} = $data[2];
}
close(SAMPLE);

#Reading Control QC
open(CONTROL, "<$controlqc") || die $!;
while(<CONTROL>){
  chomp $_;
  my @data = split("\t", $_);
  $CQC{$data[0]}{'value'} = $data[1];
  $CQC{$data[0]}{'score'} = $data[2];
}
close(CONTROL);

#Reading Overall QC
open(OVERALL, "<$overallqc") || die $!;
while(<OVERALL>){
  chomp $_;
  my @data = split("\t", $_);
  $OQC{$data[0]}{'value'} = $data[1];
  $OQC{$data[0]}{'score'} = $data[2];
}
close(OVERALL);

my $htmlheader = "<table class='results'><tr><th>DATA";
my $textheader = "DATA";
my $samplehtmlvalues = "<tr><td><center>SAMPLE</center></td>";
my $controlhtmlvalues = "<tr><td><center>CONTROL</center></td>";
my $overallhtmlvalues = "<tr><td><center>SAMPLE+CONTROL</center></td>";
my $sampletextvalues = "SAMPLE";
my $controltextvalues = "CONTROL";
my $overalltextvalues = "SAMPLE+CONTROL";

foreach (sort {$a <=> $b} keys %stats) {
  my $convertheader = $stats{$_}; $convertheader =~ s/_/ /g; #change space to underscore for txt file
  $textheader .= "\t$stats{$_}";
  $htmlheader .= "</th><th>".$convertheader;
  if (exists $SQC{$stats{$_}}) {
    $samplehtmlvalues .="<td bgcolor='".$color{$SQC{$stats{$_}}{'score'}}."'><center>".$SQC{$stats{$_}}{'value'}."</center></td>";
    $sampletextvalues .= "\t$SQC{$stats{$_}}{'value'}";
  } else {
    $samplehtmlvalues .="<td></td>"; $sampletextvalues .= "\t";
  }
  if (exists $CQC{$stats{$_}}) {
    $controlhtmlvalues .="<td bgcolor='".$color{$CQC{$stats{$_}}{'score'}}."'><center>".$CQC{$stats{$_}}{'value'}."</center></td>";
    $controltextvalues .= "\t$CQC{$stats{$_}}{'value'}";
  } else {
    $controlhtmlvalues .="<td></td>"; $controltextvalues .= "\t";
  }
  if (exists $OQC{$stats{$_}}) {
    $overallhtmlvalues .="<td bgcolor='".$color{$OQC{$stats{$_}}{'score'}}."'><center>".$OQC{$stats{$_}}{'value'}."</center></td>";
    $overalltextvalues .= "\t$OQC{$stats{$_}}{'value'}";
  } else {
    $overallhtmlvalues .="<td></td>"; $overalltextvalues .= "\t";
  }
}
$htmlheader .= "</th></tr>";
$samplehtmlvalues .= "</tr>";
$controlhtmlvalues .= "</tr>";
$overallhtmlvalues .= "</tr>";

open (OUT1, ">$htmlfile"); #creating htmlfile
print OUT1 $htmlheader, "\n", $samplehtmlvalues, "\n", $controlhtmlvalues, "\n", $overallhtmlvalues;
close (OUT1);

open (OUT2, ">$textfile"); #creating htmlfile
print OUT2 $textheader, "\n", $sampletextvalues, "\n", $controltextvalues, "\n", $overalltextvalues, "\n";
close (OUT2);

