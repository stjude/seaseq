#!/usr/bin/perl
#summary stats of all files provided

use Pod::Usage;
use strict; 
use warnings;
use File::Basename;
use Getopt::Long;

my ($help, $manual, $sppout, $bambed, $countsfile, $peaksxls, $bamflag, $rmdupflag, $bkflag, $outfile, $fastqmetrics, $fastqczip, $rosedir, $samplename);
my $usage = "perl $0 -s <spp file> -c <countsfile> -b <bambed> -px <peaks xls> -fmetric <fastq metrics> -rose <rose directory> -bamflag <bamflagstat> -fqc <fastqczipfile> [-bkflag <bklistflagstat>] [-rmdupflag <rmdupflagstat>] -outfile <outputfile>\n";
# USAGE DETAILS
#<spp file> : run_spp.R output file
#<countsfile> : number of reads overlap between mapping bam and peaks bed files
#<bambed> : mapping bam to bed file
#<peaks xls> : macs peaks xls file
#<fastq metrics> : reads fastq metrics file
#<rose directory> : rose output directory
#<fastqczipfile> : reads fastqc zip file
#<bamflagstat> : mapping bam samtools flagstat
#<bklistflagstat> : mapping bam+bklist samtools flagstat
#<rmdupflagstat> : mapping bam+rmdup samtools flagstat

GetOptions ("s|spp=s"=>\$sppout,"b|bed=s"=>\$bambed,"c|count=s"=>\$countsfile,"px|peakxls=s"=>\$peaksxls,"bamflag=s"=>\$bamflag,"rmdupflag=s"=>\$rmdupflag,"bkflag=s"=>\$bkflag,"fqc|fastqczip=s"=>\$fastqczip,"fmetric|fx=s"=>\$fastqmetrics, "rose|r=s"=>\$rosedir, "outfile|o=s"=>\$outfile);
unless ($sppout || $bambed || $countsfile || $peaksxls || $bamflag || $rmdupflag || $bkflag || $fastqczip || $fastqmetrics || $rosedir) { die $usage;}

#Filenames
my ($statsout, $htmlfile, $textfile);

#Initialize variables
my ($Uniquecnt, $Totalcnt, $Fripcnt, $FRIP, $peaks, $PBC, $NRF, $PhantomQual) = (0,0,0,0,0,0,0,0);
my (%HASH, %OvQual, $alignedpercent, $totalreads);
my $prev = "NA";

#output file name
unless ($outfile) { 
  $statsout = "summarystats-stats.out"; 
} else {
  unless ($outfile =~ /\-stats.out$/) { $statsout = fileparse($outfile, qr/(\.[\w\d]+)$/)."-stats.out"; }
  else { $statsout = $outfile; }
}
#html output file name
$htmlfile = $statsout; $htmlfile =~ s/stats.out/stats.html/; #creating html file
$textfile = $statsout; $textfile =~ s/stats.out/stats.txt/; #creating html file
open (OUT, ">$statsout"); #opening outputfile
$samplename = (split('\.', $statsout))[0];
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
} # end if fastqczip

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
  } #end if fastqczip
  
} #end if bamflag

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
if ($bambed) {
  print "Processing NRF score ...";
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
  if ($Uniquecnt >= 0 && $Totalcnt >= 0 && $Oneread >= 0) {
    $NRF = $Uniquecnt/$Totalcnt;
    $PBC = $Oneread/$Uniquecnt;
  } else { $NRF = 0; $PBC = 0; }
  print OUT "Unique Genomic Locations = $Uniquecnt\nNRF score = $NRF\nPCR Bottleneck Coefficient = $PBC\n";
  print ".. Done\n";
  
  #QCdash
  $OvQual{'Non Redundant Percent'}{'value'} = $NRF;
  $OvQual{'PCR Bottleneck'}{'value'} = $PBC;
  $OvQual{'Non Redundant Percent'}{'score'} = -2;
  $OvQual{'PCR Bottleneck'}{'score'} = -2;
  if ($OvQual{'Non Redundant Percent'}{'value'} >= 0.5) { $OvQual{'Non Redundant Percent'}{'score'} = -1; }
  if ($OvQual{'Non Redundant Percent'}{'value'} >= 0.6) { $OvQual{'Non Redundant Percent'}{'score'} = 0; }
  if ($OvQual{'Non Redundant Percent'}{'value'} >= 0.7) { $OvQual{'Non Redundant Percent'}{'score'} = 1; }
  if ($OvQual{'Non Redundant Percent'}{'value'} >= 0.8) { $OvQual{'Non Redundant Percent'}{'score'} = 2; }
  if ($OvQual{'PCR Bottleneck'}{'value'} >= 0.5) { $OvQual{'PCR Bottleneck'}{'score'} = -1; }
  if ($OvQual{'PCR Bottleneck'}{'value'} >= 0.66) { $OvQual{'PCR Bottleneck'}{'score'} = 0; }
  if ($OvQual{'PCR Bottleneck'}{'value'} >= 0.75) { $OvQual{'PCR Bottleneck'}{'score'} = 1; }
  if ($OvQual{'PCR Bottleneck'}{'value'} >= 0.9) { $OvQual{'PCR Bottleneck'}{'score'} = 2; }

}
if ($sppout) {
  #Process NSC (normalized) + RSC (relative strand cross-correlation coefficient)
  open (IN, "<$sppout"); # || die $!; 
  my ($NSC, $RSC,$PhantomQual) = (split("\t", <IN>))[8,9,10]; close(IN);
  print OUT "Normalized Strand cross-correlation Coefficient (NSC) = $NSC\n";
  print OUT "Relative Strand cross-correlation Coefficient (RSC) = $RSC\n";
  print OUT "Phantom Quality = $PhantomQual"; chop $PhantomQual;
  print ".. Done\n";
  
  #QCdash
  $OvQual{'Normalized Strand Cross-correlation'}{'value'} = $NSC;
  $OvQual{'Relative Strand Cross-correlation'}{'value'} = $RSC;
  $OvQual{'Phantom Quality'}{'value'} = $PhantomQual;
  $OvQual{'Normalized Strand Cross-correlation'}{'score'} = -2;
  $OvQual{'Relative Strand Cross-correlation'}{'score'} = -2;
  $OvQual{'Phantom Quality'}{'score'} = $PhantomQual;
  if ($OvQual{'Normalized Strand Cross-correlation'}{'value'} >= 1.1) { $OvQual{'Normalized Strand Cross-correlation'}{'score'} = -1; }
  if ($OvQual{'Normalized Strand Cross-correlation'}{'value'} >= 1.15) { $OvQual{'Normalized Strand Cross-correlation'}{'score'} = 0; }
  if ($OvQual{'Normalized Strand Cross-correlation'}{'value'} >= 1.25) { $OvQual{'Normalized Strand Cross-correlation'}{'score'} = 1; }
  if ($OvQual{'Normalized Strand Cross-correlation'}{'value'} >= 1.5) { $OvQual{'Normalized Strand Cross-correlation'}{'score'} = 2; }
  if ($OvQual{'Relative Strand Cross-correlation'}{'value'} >= 0.25) { $OvQual{'Relative Strand Cross-correlation'}{'score'} = -1; }
  if ($OvQual{'Relative Strand Cross-correlation'}{'value'} >= 0.5) { $OvQual{'Relative Strand Cross-correlation'}{'score'} = 0; }
  if ($OvQual{'Relative Strand Cross-correlation'}{'value'} >= 1) { $OvQual{'Relative Strand Cross-correlation'}{'score'} = 1; }
  if ($OvQual{'Relative Strand Cross-correlation'}{'value'} >= 1.5) { $OvQual{'Relative Strand Cross-correlation'}{'score'} = 2; }
}

if ($countsfile) {
  #FRIP score
  print "Processing FRIP score ... ";
  open(IN,"<$countsfile")|| die $!;
  while(<IN>){
    chomp;
    my @l=split("\t",$_);
    $Fripcnt+=$l[-1];
    $peaks++; 
  } close (IN);
  $FRIP = sprintf ("%.6f", ($Fripcnt/$Totalcnt));
  print OUT "Total Peaks = $peaks\nFRIP score = $FRIP\n";
  print ".. Done\n";

  #QCdash
  $OvQual{'FRiP'}{'value'} = $FRIP; $OvQual{'FRiP'}{'score'} = -2;
  if ($OvQual{'FRiP'}{'value'} >= 0.0075) {$OvQual{'FRiP'}{'score'} = -1;}
  if ($OvQual{'FRiP'}{'value'} >= 0.01) {$OvQual{'FRiP'}{'score'} = 0;}
  if ($OvQual{'FRiP'}{'value'} >= 0.02) {$OvQual{'FRiP'}{'score'} = 1;}
  if ($OvQual{'FRiP'}{'value'} >= 0.05) {$OvQual{'FRiP'}{'score'} = 2;}
  
}
  
  #Estimated fragment width & estimated tag length
if ($peaksxls) {
  print "Processing Fragment width & Tag length score ... ";
  my $fragmentwidth = `grep "\# d = " $peaksxls | head -n 1 | awk '{print \$NF}'`;
  print OUT "Estimated Fragment Width = $fragmentwidth"; chop $fragmentwidth;
  
  #Estimated Tag Length
  my $predictedtaglength = `grep "\# tag size is determined" $peaksxls | head -n 1 | awk '{print \$(NF-1)}'`;
  print OUT "Estimated Tag Length = $predictedtaglength"; chop $predictedtaglength;
  print ".. Done\n";
  #QCdash
  $OvQual{'Estimated Fragment Width'}{'value'} = $fragmentwidth; $OvQual{'Estimated Fragment Width'}{'score'} = 2;
  $OvQual{'Estimated Tag Length'}{'value'} = $predictedtaglength; $OvQual{'Estimated Tag Length'}{'score'} = 2;
  #working on the fastq metrics file to get avg read length
  if ($fastqmetrics) {
    my $avgreadlength = `tail -n 1 $fastqmetrics | awk '{print \$4}'`; chop $avgreadlength;
    if ( (($predictedtaglength - $avgreadlength) > 10) || (($predictedtaglength - $avgreadlength) < -10) ) {
      $OvQual{'Estimated Tag Length'}{'score'} = -2; #QC hash
    }
  } # end if fastqmetrics
} #end if peaksxls
 
#Stitched regions & Enhancers
if ($rosedir){
  print "Processing ROSE counts ...";
  my $enhancers = `wc -l $rosedir/unionpeaks_AllEnhancers.table.txt | awk '{print \$1-6}'`; chop $enhancers;
  my $superenhancers = `wc -l $rosedir/unionpeaks_SuperEnhancers.table.txt | awk '{print \$1-6}'`; chop $superenhancers;
  unless ($enhancers > 0 ) { $enhancers = 0; }  #making sure enhancers  is numeric
  unless ($superenhancers > 0 ) { $superenhancers = 0; } #making sure superenhancers is numeric
  print OUT "Total number of Enhancers = ", $enhancers,"\n";
  print OUT "Total number of Superenhancers = ", $superenhancers,"\n";
  print ".. Done\n";
  
  #QCdash
  $OvQual{'Enhancers'}{'value'} = $enhancers; $OvQual{'Enhancers'}{'score'} = -2;
  $OvQual{'Super Enhancers'}{'value'} = $superenhancers; $OvQual{'Super Enhancers'}{'score'} = 2;
  if ($OvQual{'Enhancers'}{'value'} > 1000) {$OvQual{'Enhancers'}{'score'} = -1;}
  if ($OvQual{'Enhancers'}{'value'} > 2000) {$OvQual{'Enhancers'}{'score'} = 0;}
  if ($OvQual{'Enhancers'}{'value'} > 5000) {$OvQual{'Enhancers'}{'score'} = 1;}
  if ($OvQual{'Enhancers'}{'value'} >= 10000) {$OvQual{'Enhancers'}{'score'} = 2;}

  if (($OvQual{'Enhancers'}{'value'} == 0) || ($OvQual{'Super Enhancers'}{'value'} == 0)) { $OvQual{'Super Enhancers'}{'score'} = -2 }
  else {
    if (($OvQual{'Super Enhancers'}{'value'}/$OvQual{'Enhancers'}{'value'}) > 0.02) {$OvQual{'Super Enhancers'}{'score'} = 1;}
    if (($OvQual{'Super Enhancers'}{'value'}/$OvQual{'Enhancers'}{'value'}) > 0.05) {$OvQual{'Super Enhancers'}{'score'} = 0;}
    if (($OvQual{'Super Enhancers'}{'value'}/$OvQual{'Enhancers'}{'value'}) > 0.1) {$OvQual{'Super Enhancers'}{'score'} = -1;}
    if (($OvQual{'Super Enhancers'}{'value'}/$OvQual{'Enhancers'}{'value'}) >= 0.2) {$OvQual{'Super Enhancers'}{'score'} = -2;}
    if ($OvQual{'Super Enhancers'}{'value'} <= 1) {$OvQual{'Super Enhancers'}{'score'} = -2;}
  }
} # end if rosedir

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

my $htmlheader = "<table border='1' cellpadding='5'><tr><th>Sample Name</th><th>"."Overall Quality";
my $textheader = "Sample Name\tOverall Quality";
my $htmlvalues = "<tr><td><center>".$samplename."</center></td>";
$htmlvalues .= "<td bgcolor='$color'><center>".$OverallQuality."</center></td>";
my $textvalues = "$samplename\t$OverallQuality";

foreach (sort keys %OvQual){
  $htmlheader .= "</th><th>".$_;
  $textheader .= "\t$_";
  $htmlvalues .="<td bgcolor='".$color{$OvQual{$_}{'score'}}."'><center>".$OvQual{$_}{'value'}."</center></td>";
  $textvalues .= "\t$OvQual{$_}{'value'}";
  print "$_\t$OvQual{$_}{'value'}\t$OvQual{$_}{'score'}\n";
}
$htmlheader .= "</th></tr>"; $htmlvalues .= "</tr>";

open (OUT2, ">$htmlfile"); #creating htmlfile
print OUT2 $htmlheader, "\n", $htmlvalues, "\n";
close (OUT2);

open (OUT3, ">$textfile"); #creating htmlfile
print OUT3 $textheader, "\n", $textvalues, "\n";
close (OUT3);
