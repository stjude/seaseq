#!/usr/bin/perl
#get region both left and right of position of interest.

use Pod::Usage;
use strict;
use warnings;
use Getopt::Long;

my ($help, $manual, $infile, $flank);
my ($verdict, $ender, $starter);

GetOptions ("bi|in=s"=>\$infile,"f|flank=i"=>\$flank);
my $usage = "perl $0 -i <bed file | gff file> -f <bp flank>\n";
unless ($infile && $flank) { die $usage; }

open(IN, $infile);
if ($infile =~ /bed$/) { $verdict = 1; $ender = 3; $starter = 0; }
elsif ($infile =~ /gff$/) { $starter = 2; $ender = 5; }

while(<IN>){
  chomp;
  my @line = split /\t/;
  if ($infile =~ /gff$/) {
    if ($line[6] =~ /\+/) {
      $verdict = 3;
    } elsif ($line[6] =~ /\-/) {
      $verdict = 4;
    }
  }
  my $end = $line[$verdict]+$flank;
  my $start = $line[$verdict]-$flank;

  print join("\t",@line[0..$starter]),"\t";
  print $start,"\t",$end,"\t";
  print join("\t",@line[$ender..$#line]),"\n";
}

