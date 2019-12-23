#!/usr/bin/perl
#summary stats of all files provided
#module needed:
#      bedtools
#      bedops
#      run_spp.R

use Pod::Usage;
use strict; 
use warnings;
use File::Basename;
use Getopt::Long;

my ($help, $manual, $inputbam, $peaksbed, $peaksxls, $bamflag, $rmdupflag, $bkflag, $outfile, $fastqmetrics, $fastqczip, $rosedir);
my $usage = "perl $0 -b <bam file> -p <peaks bed> -px <peaks xls> -fmetric <fastq metrics> -rose <rose directory> -bamflag <bamflagstat> -rmdupflag <rmdupflagstat> -bkflag <bklistflagstat> -outfile <outputfile> -fqc <fastqczipfile>\n";

GetOptions ("b|bam=s"=>\$inputbam,"p|peak=s"=>\$peaksbed,"px|peakxls=s"=>\$peaksxls,"bamflag=s"=>\$bamflag,"rmdupflag=s"=>\$rmdupflag,"bkflag=s"=>\$bkflag,"fqc|fastqczip=s"=>\$fastqczip,"fmetric|fx=s"=>\$fastqmetrics, "rose|r=s"=>\$rosedir, "outfile|o=s"=>\$outfile);
unless ($inputbam || $peaksbed || $peaksxls || $bamflag || $rmdupflag || $bkflag || $fastqczip || $fastqmetrics || $rosedir) { die $usage;}

#Filenames
my ($countsfile, $bambed, $sppout, $statsout, $tempsortbed, $htmlfile);

#Initialize variables
my ($Uniquecnt, $Totalcnt, $Fripcnt, $FRIP, $peaks, $PBC, $NRF, $PhantomQual) = (0,0,0,0,0,0,0,0);
my (%HASH, %OvQual, $alignedpercent, $totalreads);
my $prev = "NA";

#output file name
unless ($outfile) {
  if ($peaksbed) { $statsout = fileparse($peaksbed, qr/(\.[\w\d]+)?$/)."-stats.out"; } 
  else { $statsout = fileparse($inputbam, qr/(\.[\w\d]+)?$/)."-stats.out"; }
} else {
  unless ($outfile =~ /\-stats.out$/) { $statsout = fileparse($outfile, qr/(\.[\w\d]+)$/)."-stats.out"; }
  else { $statsout = $outfile; }
}
#html output file name
$htmlfile = $statsout; $htmlfile =~ s/stats.out/stats.html/; #creating html file
open (OUT, ">$statsout"); #opening outputfile

# - - - - - - - - -
# Values
# - - - - - - - - -
#working with FastQC
if ($fastqczip) {
  $totalreads = `unzip -p $fastqczip */fastqc_data.txt | grep "Total Sequences" | awk -F' ' '{print \$NF}'`;
  my $basequality = `unzip -p $fastqczip */fastqc_data.txt | grep "base sequence quality" | awk -F' ' '{print \$NF}'`;
  my $seqrep = `unzip -p $fastqczip */fastqc_data.txt | grep "Overrepresented sequences" | awk -F' ' '{print \$NF}'`;
  print OUT "Raw Reads = $totalreads"; chop $totalreads; 
  print OUT "Base Quality = $basequality"; chop $basequality;
  print OUT "Sequence Overrep = $seqrep"; chop $seqrep;
  
  #QCdash
  $OvQual{'Raw Reads'}{'value'} = $totalreads; $OvQual{'Base Quality'}{'value'} = $basequality; $OvQual{'Sequence OverRep'}{'value'} = $seqrep;
  $OvQual{'Raw Reads'}{'score'} = -2; $OvQual{'Base Quality'}{'score'} = -2; $OvQual{'Sequence OverRep'}{'score'} = -2;
  
  if ($OvQual{'Raw Reads'}{'value'} >= 15000000) { $OvQual{'Raw Reads'}{'score'} = -1; }
  if ($OvQual{'Raw Reads'}{'value'} >= 20000000) { $OvQual{'Raw Reads'}{'score'} = 0; }
  if ($OvQual{'Raw Reads'}{'value'} >= 25000000) { $OvQual{'Raw Reads'}{'score'} = 1; }
  if ($OvQual{'Raw Reads'}{'value'} >= 30000000) { $OvQual{'Raw Reads'}{'score'} = 2; }
  if ($OvQual{'Base Quality'}{'value'} eq 'warn') { $OvQual{'Base Quality'}{'score'} = 0; }
  if ($OvQual{'Base Quality'}{'value'} eq 'pass') { $OvQual{'Base Quality'}{'score'} = 2; }
  if ($OvQual{'Sequence OverRep'}{'value'} eq 'warn') { $OvQual{'Sequence OverRep'}{'score'} = 0; }
  if ($OvQual{'Sequence OverRep'}{'value'} eq 'pass') { $OvQual{'Sequence OverRep'}{'score'} = 2; }
}

#working with Flagstat
if ($bamflag) {
  my $mappedreads = `head -n 5 $bamflag | tail -n 1 | awk -F" " '{print \$1}'`;
  print OUT "Total Mapped Reads = $mappedreads"; chop $mappedreads;
  
  if ($fastqczip) {
    $alignedpercent = $mappedreads/$totalreads;
    print OUT "Aligned percentage = ",$alignedpercent, "\n";
    
    #QCdash
    $OvQual{'Aligned Percent'}{'value'} = $alignedpercent; $OvQual{'Aligned Percent'}{'score'} = -2;
    if ($OvQual{'Aligned Percent'}{'value'} >= 0.5) {$OvQual{'Aligned Percent'}{'score'} = -1;}
    if ($OvQual{'Aligned Percent'}{'value'} >= 0.6) {$OvQual{'Aligned Percent'}{'score'} = 0;}
    if ($OvQual{'Aligned Percent'}{'value'} >= 0.7) {$OvQual{'Aligned Percent'}{'score'} = 1;}
    if ($OvQual{'Aligned Percent'}{'value'} >= 0.8) {$OvQual{'Aligned Percent'}{'score'} = 2;}
  }
}

if ($bkflag) {
  my $bklistreads = `head -n 5 $bkflag | tail -n 1 | awk -F" " '{print \$1}'`;
  print OUT "Total After RM BlackList regions = ",$bklistreads;
}

if ($rmdupflag) {
  my $rmdupreads = `head -n 5 $rmdupflag | tail -n 1 | awk -F" " '{print \$1}'`;
  print OUT "Total After RM Duplicates = ",$rmdupreads;
}


#PROCESS NRF (Non-Redundant Fraction)
#the ratio between the number of positions in the genome that uniquely mappable reads map to and the total number of uniquely mappable reads
if ($inputbam) {
  print "Processing NRF score ...";
  unless ($inputbam =~ /\.bam$/) { die $!; }
  $bambed = fileparse($inputbam, qr/(\.bam)?$/)."-bam2bed.bed";
  $sppout = fileparse($inputbam, qr/(\.bam)?$/)."-spp.out";

  `bamToBed -i $inputbam > $bambed`;
  open(IN, "<$bambed") || die $!;
  while(<IN>){
    chomp;
    my @line = split("\t",$_); 
    $Totalcnt++; #from bam file mapped reads
    my $t = join("_",@line[0..2]);
    #$Uniquecnt++ unless($t eq $prev);
    unless ($t eq $prev) {
      $Uniquecnt++;
    } else {
      $HASH{$t} = 1;
    }
    $prev=$t;
  } close (IN);
 
  my $Oneread = $Uniquecnt - scalar keys %HASH;
  $NRF = $Uniquecnt/$Totalcnt;
  $PBC = $Oneread/$Uniquecnt;
  print OUT "Unique Genomic Locations = $Uniquecnt\nNRF score = $NRF\nPCR Bottleneck Coefficient = $PBC\n";
  print ".. Done\n";


  #Process NSC (normalized) + RSC (relative strand cross-correlation coefficient)
  print "\nrun_spp.R -c=$inputbam -savp -out=$sppout 1>outfile.fake 2>outfile.fake;\n";
  print "Processing NSC and RSC ...";
  `run_spp.R -c=$inputbam -savp -out=$sppout 1>outfile.fake 2>outfile.fake;`;
  open (IN, "<$sppout"); # || die $!; 
  my ($NSC, $RSC,$PhantomQual) = (split("\t", <IN>))[8,9,10]; close(IN);
  print OUT "Normalized Strand cross-correlation Coefficient (NSC) = $NSC\n";
  print OUT "Relative Strand cross-correlation Coefficient (RSC) = $RSC\n";
  print OUT "Phantom Quality = $PhantomQual"; chop $PhantomQual;
  print ".. Done\n";
  
  #QCdash
  $OvQual{'Non Redundant Percent'}{'value'} = $NRF;
  $OvQual{'PCR Bottleneck'}{'value'} = $PBC;
  $OvQual{'Normalized Strand Cross-correlation'}{'value'} = $NSC;
  $OvQual{'Relative Strand Cross-correlation'}{'value'} = $RSC;
  $OvQual{'Phantom Quality'}{'value'} = $PhantomQual;
  $OvQual{'Non Redundant Percent'}{'score'} = -2;
  $OvQual{'PCR Bottleneck'}{'score'} = -2;
  $OvQual{'Normalized Strand Cross-correlation'}{'score'} = -2;
  $OvQual{'Relative Strand Cross-correlation'}{'score'} = -2;
  $OvQual{'Phantom Quality'}{'score'} = $PhantomQual;
  if ($OvQual{'Non Redundant Percent'}{'value'} >= 0.5) { $OvQual{'Non Redundant Percent'}{'score'} = -1; }
  if ($OvQual{'Non Redundant Percent'}{'value'} >= 0.6) { $OvQual{'Non Redundant Percent'}{'score'} = 0; }
  if ($OvQual{'Non Redundant Percent'}{'value'} >= 0.7) { $OvQual{'Non Redundant Percent'}{'score'} = 1; }
  if ($OvQual{'Non Redundant Percent'}{'value'} >= 0.8) { $OvQual{'Non Redundant Percent'}{'score'} = 2; }
  if ($OvQual{'PCR Bottleneck'}{'value'} >= 0.5) { $OvQual{'PCR Bottleneck'}{'score'} = -1; }
  if ($OvQual{'PCR Bottleneck'}{'value'} >= 0.66) { $OvQual{'PCR Bottleneck'}{'score'} = 0; }
  if ($OvQual{'PCR Bottleneck'}{'value'} >= 0.75) { $OvQual{'PCR Bottleneck'}{'score'} = 1; }
  if ($OvQual{'PCR Bottleneck'}{'value'} >= 0.9) { $OvQual{'PCR Bottleneck'}{'score'} = 2; }
  if ($OvQual{'Normalized Strand Cross-correlation'}{'value'} >= 1.1) { $OvQual{'Normalized Strand Cross-correlation'}{'score'} = -1; }
  if ($OvQual{'Normalized Strand Cross-correlation'}{'value'} >= 1.15) { $OvQual{'Normalized Strand Cross-correlation'}{'score'} = 0; }
  if ($OvQual{'Normalized Strand Cross-correlation'}{'value'} >= 1.25) { $OvQual{'Normalized Strand Cross-correlation'}{'score'} = 1; }
  if ($OvQual{'Normalized Strand Cross-correlation'}{'value'} >= 1.5) { $OvQual{'Normalized Strand Cross-correlation'}{'score'} = 2; }
  if ($OvQual{'Relative Strand Cross-correlation'}{'value'} >= 0.25) { $OvQual{'Relative Strand Cross-correlation'}{'score'} = -1; }
  if ($OvQual{'Relative Strand Cross-correlation'}{'value'} >= 0.5) { $OvQual{'Relative Strand Cross-correlation'}{'score'} = 0; }
  if ($OvQual{'Relative Strand Cross-correlation'}{'value'} >= 1) { $OvQual{'Relative Strand Cross-correlation'}{'score'} = 1; }
  if ($OvQual{'Relative Strand Cross-correlation'}{'value'} >= 1.5) { $OvQual{'Relative Strand Cross-correlation'}{'score'} = 2; }
}

if ($peaksbed && $inputbam) {
#FRIP score
print "Processing FRIP score ... ";
  $countsfile = fileparse($peaksbed, qr/(\.bed)?$/)."-out.txt"; 
  $tempsortbed = fileparse($bambed, qr/(\.bed)?$/).".sorted.bed";
  `sort-bed $bambed > $tempsortbed`;
  `intersectBed -sorted -a $peaksbed -b $tempsortbed -c > $countsfile`;

  open(IN,"<$countsfile")|| die $!;
  while(<IN>){
    chomp;
    my @l=split("\t",$_);
    $Fripcnt+=$l[-1];
    $peaks++; 
  } close (IN);
  $FRIP = sprintf ("%.6f", ($Fripcnt/$Totalcnt));
  print OUT "Total Peaks = $peaks\nFRIP score = $FRIP\n";

  #QCdash
  $OvQual{'FRiP'}{'value'} = $FRIP; $OvQual{'FRiP'}{'score'} = -2;
  if ($OvQual{'FRiP'}{'value'} >= 0.0075) {$OvQual{'FRiP'}{'score'} = -1;}
  if ($OvQual{'FRiP'}{'value'} >= 0.01) {$OvQual{'FRiP'}{'score'} = 0;}
  if ($OvQual{'FRiP'}{'value'} >= 0.02) {$OvQual{'FRiP'}{'score'} = 1;}
  if ($OvQual{'FRiP'}{'value'} >= 0.05) {$OvQual{'FRiP'}{'score'} = 2;}  
  
  #Estimated fragment width & estimated tag length
  if ($peaksxls) {
    my $fragmentwidth = `grep "\# d = " $peaksxls | head -n 1 | awk '{print \$NF}'`;
    print OUT "Estimated Fragment Width = $fragmentwidth"; chop $fragmentwidth;
    
    #Estimated Tag Length
    my $predictedtaglength = `grep "\# tag size is determined" $peaksxls | head -n 1 | awk '{print \$(NF-1)}'`;
    print OUT "Estimated Tag Length = $predictedtaglength"; chop $predictedtaglength;
    
    #QCdash
    $OvQual{'Estimated Fragment Width'}{'value'} = $fragmentwidth; $OvQual{'Estimated Fragment Width'}{'score'} = 2;
    $OvQual{'Estimated Tag Length'}{'value'} = $predictedtaglength; $OvQual{'Estimated Tag Length'}{'score'} = 2;
    #working on the fastq metrics file to get avg read length
    if ($fastqmetrics) {
      my $avgreadlength = `tail -n 1 $fastqmetrics | awk '{print \$4}'`; chop $avgreadlength;
      if ( (($predictedtaglength - $avgreadlength) > 10) || (($predictedtaglength - $avgreadlength) < -10) ) { $OvQual{'Estimated Tag Length'}{'score'} = -2; }
    }
  }
  print ".. Done\n";
}

if ($rosedir){
  print "Processing ROSE counts ...";
  my $enhancers = `wc -l $rosedir/unionpeaks_AllEnhancers.table.txt | awk '{print \$1-6}'`; chop $enhancers;
  my $superenhancers = `wc -l $rosedir/unionpeaks_SuperEnhancers.table.txt | awk '{print \$1-6}'`; chop $superenhancers;
  print OUT "Total number of Enhancers = ", $enhancers;
  print OUT "Total number of Superenhancers = ", $superenhancers;
  print ".. Done\n";
  
  #QCdash
  $OvQual{'Enhancers'}{'value'} = $enhancers; $OvQual{'Enhancers'}{'score'} = -2;
  $OvQual{'Super Enhancers'}{'value'} = $superenhancers; $OvQual{'Super Enhancers'}{'score'} = 2;
  if ($OvQual{'Enhancers'}{'value'} > 1000) {$OvQual{'Enhancers'}{'score'} = -1;}
  if ($OvQual{'Enhancers'}{'value'} > 2000) {$OvQual{'Enhancers'}{'score'} = 0;}
  if ($OvQual{'Enhancers'}{'value'} > 5000) {$OvQual{'Enhancers'}{'score'} = 1;}
  if ($OvQual{'Enhancers'}{'value'} >= 10000) {$OvQual{'Enhancers'}{'score'} = 2;}
  if (($OvQual{'Super Enhancers'}{'value'}/$OvQual{'Enhancers'}{'value'}) > 0.02) {$OvQual{'Super Enhancers'}{'score'} = 1;}
  if (($OvQual{'Super Enhancers'}{'value'}/$OvQual{'Enhancers'}{'value'}) > 0.05) {$OvQual{'Super Enhancers'}{'score'} = 0;}
  if (($OvQual{'Super Enhancers'}{'value'}/$OvQual{'Enhancers'}{'value'}) > 0.1) {$OvQual{'Super Enhancers'}{'score'} = -1;}
  if (($OvQual{'Super Enhancers'}{'value'}/$OvQual{'Enhancers'}{'value'}) >= 0.2) {$OvQual{'Super Enhancers'}{'score'} = -2;}
  if ($OvQual{'Super Enhancers'}{'value'} <= 1) {$OvQual{'Super Enhancers'}{'score'} = -2;}
}

close(OUT);

#Processing QC html dashboard
my ($count,$totalscore) = (0,0); 
foreach (keys %OvQual){
  $count++;
  $totalscore += $OvQual{$_}{'score'};
}
#color names
my %color = ( "-2" => "#FF0000", "-1" => "#FF8C00", "0" => "#FFFF00", "1" => "#ADFF2F", "2" => "#008000" ); #red #orangered #yellow #greenyellow #green

my $OverallQuality = sprintf ("%.4f", ($totalscore/$count)); my $color = $color{(sprintf ("%.0f", ($totalscore/$count)))};

my $header = "<table border='1' cellpadding='5'><tr><th>"."Overall Quality";
my $values = "<tr><td bgcolor='$color'><center>".$OverallQuality."</center></td>";
foreach (sort keys %OvQual){
  $header .= "</th><th>".$_;
  $values .="<td bgcolor='".$color{$OvQual{$_}{'score'}}."'><center>".$OvQual{$_}{'value'}."</center></td>";
  print "$_\t$OvQual{$_}{'value'}\t$OvQual{$_}{'score'}\n";
}
$header .= "</th></tr>"; $values .= "</tr>";
my $end = "</table>";

open (OUT2, ">$htmlfile"); #creating htmlfile
print OUT2 $header, "\n", $values, "\n", $end;
close (OUT2);
