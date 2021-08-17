#!/usr/bin/perl
#summary stats of all files provided

use Pod::Usage;
use strict; 
use warnings;
use File::Basename;
use Getopt::Long;

my ($help, $manual, $sppout, $bambed, $countsfile, $peaksxls, $bamflag, $rmdupflag,
    $bkflag, $outfile, $fastqmetrics, $fastqczip, $rose_stitched, $rose_superstitched,
    $samplename, $cbamflag, $crmdupflag, $cbkflag, $cfastqczip, $ccountsfile, $cbambed,
    $cpeaksxls);
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

GetOptions (
            "s|spp=s"=>\$sppout, "b|bed=s"=>\$bambed, "c|count=s"=>\$countsfile,
            "px|peakxls=s"=>\$peaksxls, "bamflag=s"=>\$bamflag, "rmdupflag=s"=>\$rmdupflag,
            "bkflag=s"=>\$bkflag, "fqc|fastqczip=s"=>\$fastqczip, "fmetric|fx=s"=>\$fastqmetrics,
            "roseenhancers|re=s"=>\$rose_stitched, "rosesuper|rs=s"=>\$rose_superstitched,
            "cbamflag=s"=>\$cbamflag, "crmdupflag=s"=>\$crmdupflag, "cbkflag=s"=>\$cbkflag,
            "cfqc|cfastqczip=s"=>\$cfastqczip, "ccount=s"=>\$ccountsfile, "cpx|cpeaksxls=s"=>\$cpeaksxls,
            "cb|cbed=s"=>\$cbambed, "outfile|o=s"=>\$outfile);
unless ($sppout || $bambed || $countsfile || $peaksxls || $bamflag || $rmdupflag || $bkflag || $fastqczip || $fastqmetrics || $rose_stitched || $rose_superstitched) { die $usage;}

#Filenames
my ($statsout, $htmlfile, $textfile, $configfile);

#Initialize variables
my ($Uniquecnt, $Totalcnt, $Fripcnt, $FRIP, $peaks, $PBC, $NRF, $PhantomQual) = (0,0,0,0,0,0,0,0);
my ($cUniquecnt, $cTotalcnt, $cFripcnt, $cFRIP, $cpeaks) = (0,0,0,0,0);
my (%HASH, %cHASH, %OVAL, %CONTROL_OVAL, %OvQual, %Control_OvQual, $alignedpercent, $totalreads, $calignedpercent, $ctotalreads);
my $prev = "NA";
my $output_counter = 7;

#output file name
unless ($outfile) { 
  $statsout = "summarystats-stats.csv"; 
} else {
  unless ($outfile =~ /\-stats.csv$/) { $statsout = fileparse($outfile, qr/(\.[\w\d]+)$/)."-stats.csv"; }
  else { $statsout = $outfile; }
}
#html output file name
$htmlfile = $statsout; $htmlfile =~ s/stats.csv/stats.htmlx/; #creating html file
$textfile = $statsout; $textfile =~ s/stats.csv/stats.txt/; #creating html file
$configfile = $statsout; $configfile =~ s/stats.csv/config.ml/; #creating config file
open (OUT, ">$statsout"); #opening outputfile
#$samplename = ((split('-stats', $statsout))[0]);
$samplename = (split('.sorted', ((split('-stats', $statsout))[0])))[0];
# - - - - - - - - -
# Values
# - - - - - - - - -
#working with FastQC
if ($fastqczip) {
  $totalreads = `unzip -p $fastqczip */fastqc_data.txt | grep "Total Sequences" | awk -F' ' '{print \$NF}'`;
  my $basequality = `unzip -p $fastqczip */fastqc_data.txt | grep "base sequence quality" | awk -F' ' '{print \$NF}'`;
  my $seqrep = `unzip -p $fastqczip */fastqc_data.txt | grep "Overrepresented sequences" | awk -F' ' '{print \$NF}'`;
  print OUT "Raw Reads,$totalreads"; chop $totalreads;
  print OUT "Base Quality,$basequality"; chop $basequality;
  print OUT "Sequence Diversity,$seqrep"; chop $seqrep;
  
  #QCdash
  $OVAL{0} = 'Raw Reads';
  $OVAL{2} = 'Base Quality';
  $OVAL{3} = 'Sequence Diversity';
  
  $OvQual{'Raw Reads'}{'value'} = $totalreads; $OvQual{'Base Quality'}{'value'} = $basequality; $OvQual{'Sequence Diversity'}{'value'} = $seqrep;
  $OvQual{'Raw Reads'}{'score'} = -2; $OvQual{'Base Quality'}{'score'} = -2; $OvQual{'Sequence Diversity'}{'score'} = -2;
  
  if ($OvQual{'Raw Reads'}{'value'} >= 15000000) { $OvQual{'Raw Reads'}{'score'} = -1; }
  if ($OvQual{'Raw Reads'}{'value'} >= 20000000) { $OvQual{'Raw Reads'}{'score'} = 0; }
  if ($OvQual{'Raw Reads'}{'value'} >= 25000000) { $OvQual{'Raw Reads'}{'score'} = 1; }
  if ($OvQual{'Raw Reads'}{'value'} >= 30000000) { $OvQual{'Raw Reads'}{'score'} = 2; }
  if ($OvQual{'Base Quality'}{'value'} eq 'warn') { $OvQual{'Base Quality'}{'score'} = 0; }
  if ($OvQual{'Base Quality'}{'value'} eq 'pass') { $OvQual{'Base Quality'}{'score'} = 2; }
  if ($OvQual{'Sequence Diversity'}{'value'} eq 'warn') { $OvQual{'Sequence Diversity'}{'score'} = 0; }
  if ($OvQual{'Sequence Diversity'}{'value'} eq 'pass') { $OvQual{'Sequence Diversity'}{'score'} = 2; }
} # end if fastqczip

#ReadLength
unless ($peaksxls) {
  if ($fastqmetrics) {
    $OVAL{1} = 'Read Length';
    my $avgreadlength = `tail -n 1 $fastqmetrics | awk '{print \$4}'`; chop $avgreadlength;
    $OvQual{'Read Length'}{'value'} = $peaks; $OvQual{'Read Length'}{'score'} = 2;
    print OUT "Read Length,$avgreadlength\n";
  } # end if fastqmetrics
}

#working with Flagstat
if ($bamflag) {
  my $mappedreads = `head -n 5 $bamflag | tail -n 1 | awk -F" " '{print \$1}'`;
  print OUT "Total Mapped Reads,$mappedreads"; chop $mappedreads;
  
  if ($fastqczip) {
    $alignedpercent = sprintf ("%.3f", ($mappedreads/$totalreads * 100));
    print OUT "Aligned percentage,",$alignedpercent, "\n";
    
    #QCdash
    $OVAL{4} = 'Aligned Percent';
    $OvQual{'Aligned Percent'}{'value'} = $alignedpercent; $OvQual{'Aligned Percent'}{'score'} = -2;
    if ($OvQual{'Aligned Percent'}{'value'} >= 50) {$OvQual{'Aligned Percent'}{'score'} = -1;}
    if ($OvQual{'Aligned Percent'}{'value'} >= 60) {$OvQual{'Aligned Percent'}{'score'} = 0;}
    if ($OvQual{'Aligned Percent'}{'value'} >= 70) {$OvQual{'Aligned Percent'}{'score'} = 1;}
    if ($OvQual{'Aligned Percent'}{'value'} >= 80) {$OvQual{'Aligned Percent'}{'score'} = 2;}
  } #end if fastqczip
  
} #end if bamflag

if ($bkflag) {
  my $bklistreads = `head -n 5 $bkflag | tail -n 1 | awk -F" " '{print \$1}'`;
  print OUT "Total After RM BlackList regions,",$bklistreads;
}

if ($rmdupflag) {
  my $rmdupreads = `head -n 5 $rmdupflag | tail -n 1 | awk -F" " '{print \$1}'`;
  print OUT "Total After RM Duplicates,",$rmdupreads;
}

#Working with control
#working with FastQC
if ($cfastqczip) {
  $ctotalreads = `unzip -p $cfastqczip */fastqc_data.txt | grep "Total Sequences" | awk -F' ' '{print \$NF}'`;
  my $cbasequality = `unzip -p $fastqczip */fastqc_data.txt | grep "base sequence quality" | awk -F' ' '{print \$NF}'`;
  my $cseqrep = `unzip -p $cfastqczip */fastqc_data.txt | grep "Overrepresented sequences" | awk -F' ' '{print \$NF}'`;
  print OUT "Control BAM : Raw Reads,$ctotalreads"; chop $ctotalreads;
  print OUT "Control BAM : Base Quality,$cbasequality"; chop $cbasequality;
  print OUT "Control BAM : Sequence Diversity,$cseqrep"; chop $cseqrep;
  
  #QCdash
  $CONTROL_OVAL{0} = 'Raw Reads';
  $CONTROL_OVAL{2} = 'Base Quality';
  $CONTROL_OVAL{3} = 'Sequence Diversity';
  
  $Control_OvQual{'Raw Reads'}{'value'} = $ctotalreads; $Control_OvQual{'Base Quality'}{'value'} = $cbasequality; $Control_OvQual{'Sequence Diversity'}{'value'} = $cseqrep;
  $Control_OvQual{'Raw Reads'}{'score'} = -2; $Control_OvQual{'Base Quality'}{'score'} = -2; $Control_OvQual{'Sequence Diversity'}{'score'} = -2;
  
  if ($Control_OvQual{'Raw Reads'}{'value'} >= 15000000) { $Control_OvQual{'Raw Reads'}{'score'} = -1; }
  if ($Control_OvQual{'Raw Reads'}{'value'} >= 20000000) { $Control_OvQual{'Raw Reads'}{'score'} = 0; }
  if ($Control_OvQual{'Raw Reads'}{'value'} >= 25000000) { $Control_OvQual{'Raw Reads'}{'score'} = 1; }
  if ($Control_OvQual{'Raw Reads'}{'value'} >= 30000000) { $Control_OvQual{'Raw Reads'}{'score'} = 2; }
  if ($Control_OvQual{'Base Quality'}{'value'} eq 'warn') { $Control_OvQual{'Base Quality'}{'score'} = 0; }
  if ($Control_OvQual{'Base Quality'}{'value'} eq 'pass') { $Control_OvQual{'Base Quality'}{'score'} = 2; }
  if ($Control_OvQual{'Sequence Diversity'}{'value'} eq 'warn') { $Control_OvQual{'Sequence Diversity'}{'score'} = 0; }
  if ($Control_OvQual{'Sequence Diversity'}{'value'} eq 'pass') { $Control_OvQual{'Sequence Diversity'}{'score'} = 2; }
} # end if fastqczip

#Working with Flagstat
if ($cbamflag) {
  my $cmappedreads = `head -n 5 $cbamflag | tail -n 1 | awk -F" " '{print \$1}'`;
  print OUT "Control BAM : Total Mapped Reads,$cmappedreads"; chop $cmappedreads;
  
  if ($cfastqczip) {
    $calignedpercent = sprintf ("%.3f", ($cmappedreads/$ctotalreads * 100));
    print OUT "Control BAM : Aligned percentage,",$calignedpercent, "\n";
    
    #QCdash
    $CONTROL_OVAL{4} = 'Aligned Percent';
    $Control_OvQual{'Aligned Percent'}{'value'} = $calignedpercent; $Control_OvQual{'Aligned Percent'}{'score'} = -2;
    if ($Control_OvQual{'Aligned Percent'}{'value'} >= 50) {$Control_OvQual{'Aligned Percent'}{'score'} = -1;}
    if ($Control_OvQual{'Aligned Percent'}{'value'} >= 60) {$Control_OvQual{'Aligned Percent'}{'score'} = 0;}
    if ($Control_OvQual{'Aligned Percent'}{'value'} >= 70) {$Control_OvQual{'Aligned Percent'}{'score'} = 1;}
    if ($Control_OvQual{'Aligned Percent'}{'value'} >= 80) {$Control_OvQual{'Aligned Percent'}{'score'} = 2;}
  } #end if cfastqczip
  
} #end if bamflag

if ($cbkflag) {
  my $cbklistreads = `head -n 5 $cbkflag | tail -n 1 | awk -F" " '{print \$1}'`;
  print OUT "Control BAM : Total After RM BlackList regions,",$cbklistreads;
}

if ($crmdupflag) {
  my $crmdupreads = `head -n 5 $crmdupflag | tail -n 1 | awk -F" " '{print \$1}'`;
  print OUT "Control BAM : Total After RM Duplicates,",$crmdupreads;
}
#end of working with control

#PROCESS NRF (Non-Redundant Fraction)
#the ratio between the number of positions in the genome that uniquely mappable reads map to and the total number of uniquely mappable reads
if ($bambed) {
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
    $NRF = sprintf ("%.4f", ($Uniquecnt/$Totalcnt));
    $PBC = sprintf ("%.4f", ($Oneread/$Uniquecnt));
  } else { $NRF = 0; $PBC = 0; }
  print OUT "Unique Genomic Locations,$Uniquecnt\nNRF score,$NRF\nPCR Bottleneck Coefficient,$PBC\n";
  
  #QCdash
  $OVAL{$output_counter++} = 'NRF';
  $OVAL{$output_counter++} = 'PBC';

  $OvQual{'NRF'}{'value'} = $NRF;
  $OvQual{'PBC'}{'value'} = $PBC;
  $OvQual{'NRF'}{'score'} = -2;
  $OvQual{'PBC'}{'score'} = -2;
  if ($OvQual{'NRF'}{'value'} >= 0.5) { $OvQual{'NRF'}{'score'} = -1; }
  if ($OvQual{'NRF'}{'value'} >= 0.6) { $OvQual{'NRF'}{'score'} = 0; }
  if ($OvQual{'NRF'}{'value'} >= 0.7) { $OvQual{'NRF'}{'score'} = 1; }
  if ($OvQual{'NRF'}{'value'} >= 0.8) { $OvQual{'NRF'}{'score'} = 2; }
  if ($OvQual{'PBC'}{'value'} >= 0.5) { $OvQual{'PBC'}{'score'} = -1; }
  if ($OvQual{'PBC'}{'value'} >= 0.66) { $OvQual{'PBC'}{'score'} = 0; }
  if ($OvQual{'PBC'}{'value'} >= 0.75) { $OvQual{'PBC'}{'score'} = 1; }
  if ($OvQual{'PBC'}{'value'} >= 0.9) { $OvQual{'PBC'}{'score'} = 2; }
}

if ($sppout) {
  #Process NSC (normalized) + RSC (relative strand cross-correlation coefficient)
  open (IN, "<$sppout"); # || die $!; 
  my ($NSC, $RSC,$PhantomQual) = (split("\t", <IN>))[8,9,10]; close(IN);
  $NSC = sprintf ("%.4f", $NSC);
  $RSC = sprintf ("%.4f", $RSC);
  print OUT "Normalized Strand cross-correlation Coefficient (NSC),$NSC\n";
  print OUT "Relative Strand cross-correlation Coefficient (RSC),$RSC\n";
  print OUT "Phantom Quality,$PhantomQual"; chop $PhantomQual;
  
  #QCdash
  $OVAL{$output_counter++} = 'NSC';
  $OVAL{$output_counter++} = 'RSC';

  $OvQual{'NSC'}{'value'} = $NSC;
  $OvQual{'RSC'}{'value'} = $RSC;
  #$OvQual{'Phantom Quality'}{'value'} = $PhantomQual; #removed from Overall Quality because it's scale of NSC & RSC (redundant)
  $OvQual{'NSC'}{'score'} = -2;
  $OvQual{'RSC'}{'score'} = -2;
  #$OvQual{'Phantom Quality'}{'score'} = $PhantomQual; #removed from Overall Quality
  if ($OvQual{'NSC'}{'value'} >= 1.045) { $OvQual{'NSC'}{'score'} = 2; }
  if ($OvQual{'RSC'}{'value'} >= 0.75) { $OvQual{'RSC'}{'score'} = 1; }
  if ($OvQual{'RSC'}{'value'} >= 1) { $OvQual{'RSC'}{'score'} = 2; }
}

if ($countsfile) {
  #FRIP score
  print "Processing FRIP score ... \n";
  open(IN,"<$countsfile")|| die $!;
  while(<IN>){
    chomp;
    my @l=split("\t",$_);
    $Fripcnt+=$l[-1];
    $peaks++; 
  } close (IN);
  $FRIP = sprintf ("%.4f", ($Fripcnt/$Totalcnt));
  print OUT "Total Peaks,$peaks\nFRIP score,$FRIP\n";

  #QCdash
  $OVAL{$output_counter++} = 'FRiP';
  $OVAL{$output_counter++} = 'Total Peaks';

  $OvQual{'Total Peaks'}{'value'} = $peaks; $OvQual{'Total Peaks'}{'score'} = -2;
  if ($OvQual{'Total Peaks'}{'value'} > 1000) {$OvQual{'Total Peaks'}{'score'} = -1;}
  if ($OvQual{'Total Peaks'}{'value'} > 2000) {$OvQual{'Total Peaks'}{'score'} = 0;}
  if ($OvQual{'Total Peaks'}{'value'} > 5000) {$OvQual{'Total Peaks'}{'score'} = 1;}
  if ($OvQual{'Total Peaks'}{'value'} >= 10000) {$OvQual{'Total Peaks'}{'score'} = 2;}
  $OvQual{'FRiP'}{'value'} = $FRIP; $OvQual{'FRiP'}{'score'} = -2;
  if ($OvQual{'FRiP'}{'value'} >= 0.0075) {$OvQual{'FRiP'}{'score'} = -1;}
  if ($OvQual{'FRiP'}{'value'} >= 0.01) {$OvQual{'FRiP'}{'score'} = 0;}
  if ($OvQual{'FRiP'}{'value'} >= 0.02) {$OvQual{'FRiP'}{'score'} = 1;}
  if ($OvQual{'FRiP'}{'value'} >= 0.05) {$OvQual{'FRiP'}{'score'} = 2;}
  
}

#working with Control Peaks
if ($ccountsfile) {
  #x3x
  
  if ($cbambed) {
    open(IN, "<$cbambed") || die $!;
    while(<IN>){
      chomp;
      my @line = split("\t",$_); 
      $cTotalcnt++; #from bam file mapped reads
      my $t = join("_",@line[0..2]);
      #$Uniquecnt++ unless($t eq $prev);
      unless ($t eq $prev) {
        $cUniquecnt++;
      } else {
        $cHASH{$t} = 1;
      }
      $prev=$t;
    } close (IN);
   
    my $cOneread = $cUniquecnt - scalar keys %cHASH;
    if ($cUniquecnt >= 0 && $cTotalcnt >= 0 && $cOneread >= 0) {
      $NRF = sprintf ("%.4f", ($cUniquecnt/$cTotalcnt));
      $PBC = sprintf ("%.4f", ($cOneread/$cUniquecnt));
    } else { $NRF = 0; $PBC = 0; }
    print OUT "Unique Genomic Locations,$cUniquecnt\nNRF score,$NRF\nPCR Bottleneck Coefficient,$PBC\n";
    
    #QCdash
    $CONTROL_OVAL{$output_counter++} = 'NRF';
    $CONTROL_OVAL{$output_counter++} = 'PBC';
  
    $Control_OvQual{'NRF'}{'value'} = $NRF;
    $Control_OvQual{'PBC'}{'value'} = $PBC;
    $Control_OvQual{'NRF'}{'score'} = -2;
    $Control_OvQual{'PBC'}{'score'} = -2;
    if ($Control_OvQual{'NRF'}{'value'} >= 0.5) { $Control_OvQual{'NRF'}{'score'} = -1; }
    if ($Control_OvQual{'NRF'}{'value'} >= 0.6) { $Control_OvQual{'NRF'}{'score'} = 0; }
    if ($Control_OvQual{'NRF'}{'value'} >= 0.7) { $Control_OvQual{'NRF'}{'score'} = 1; }
    if ($Control_OvQual{'NRF'}{'value'} >= 0.8) { $Control_OvQual{'NRF'}{'score'} = 2; }
    if ($Control_OvQual{'PBC'}{'value'} >= 0.5) { $Control_OvQual{'PBC'}{'score'} = -1; }
    if ($Control_OvQual{'PBC'}{'value'} >= 0.66) { $Control_OvQual{'PBC'}{'score'} = 0; }
    if ($Control_OvQual{'PBC'}{'value'} >= 0.75) { $Control_OvQual{'PBC'}{'score'} = 1; }
    if ($Control_OvQual{'PBC'}{'value'} >= 0.9) { $Control_OvQual{'PBC'}{'score'} = 2; }
  }
    
  #FRIP score
  print "Processing FRIP score for control ... \n";
  open(IN,"<$ccountsfile")|| die $!;
  while(<IN>){
    chomp;
    my @l=split("\t",$_);
    $cFripcnt+=$l[-1];
    $cpeaks++; 
  } close (IN);
  $cFRIP = sprintf ("%.4f", ($cFripcnt/$cTotalcnt));
  print OUT "Control: FRIP score,$cFRIP\n";

  #QCdash
  $CONTROL_OVAL{$output_counter++} = 'FRiP';

  $Control_OvQual{'FRiP'}{'value'} = $cFRIP; $Control_OvQual{'FRiP'}{'score'} = -2;
  if ($Control_OvQual{'FRiP'}{'value'} >= 0.0075) {$Control_OvQual{'FRiP'}{'score'} = -1;}
  if ($Control_OvQual{'FRiP'}{'value'} >= 0.01) {$Control_OvQual{'FRiP'}{'score'} = 0;}
  if ($Control_OvQual{'FRiP'}{'value'} >= 0.02) {$Control_OvQual{'FRiP'}{'score'} = 1;}
  if ($Control_OvQual{'FRiP'}{'value'} >= 0.05) {$Control_OvQual{'FRiP'}{'score'} = 2;}
  #x3X
  
  #QCdash
  print OUT "Control: Total Peaks,$cpeaks\n";
  $CONTROL_OVAL{$output_counter++} = 'Total Peaks';

  $Control_OvQual{'Total Peaks'}{'value'} = $cpeaks; $Control_OvQual{'Total Peaks'}{'score'} = -2;
  if ($Control_OvQual{'Total Peaks'}{'value'} > 1000) {$Control_OvQual{'Total Peaks'}{'score'} = -1;}
  if ($Control_OvQual{'Total Peaks'}{'value'} > 2000) {$Control_OvQual{'Total Peaks'}{'score'} = 0;}
  if ($Control_OvQual{'Total Peaks'}{'value'} > 5000) {$Control_OvQual{'Total Peaks'}{'score'} = 1;}
  if ($Control_OvQual{'Total Peaks'}{'value'} >= 10000) {$Control_OvQual{'Total Peaks'}{'score'} = 2;}
  
}

#Estimated fragment width & estimated tag length
if ($peaksxls) {
  print "Processing Fragment width & Tag length score ... \n";
  my $fragmentwidth = `grep "\# d = " $peaksxls | head -n 1 | awk '{print \$NF}'`;
  print OUT "Estimated Fragment Width,$fragmentwidth"; chop $fragmentwidth;
  
  #Estimated Tag Length
  my $predictedtaglength = `grep "\# tag size is determined" $peaksxls | head -n 1 | awk '{print \$(NF-1)}'`;
  print OUT "Estimated Tag Length,$predictedtaglength"; chop $predictedtaglength;

  #QCdash
  $OVAL{5} = 'Estimated Fragment Width';
  $OVAL{6} = 'Estimated Tag Length';

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

#Control Estimated fragment width & estimated tag length
if ($cpeaksxls) {
  print "Processing Fragment width & Tag length score for control ... \n";
  my $cfragmentwidth = `grep "\# d = " $cpeaksxls | head -n 1 | awk '{print \$NF}'`;
  print OUT "Control: Estimated Fragment Width,$cfragmentwidth"; chop $cfragmentwidth;
  
  #Estimated Tag Length
  my $cpredictedtaglength = `grep "\# tag size is determined" $cpeaksxls | head -n 1 | awk '{print \$(NF-1)}'`;
  print OUT "Control: Estimated Tag Length,$cpredictedtaglength"; chop $cpredictedtaglength;

  #QCdash
  $CONTROL_OVAL{5} = 'Estimated Fragment Width';
  $CONTROL_OVAL{6} = 'Estimated Tag Length';

  $Control_OvQual{'Estimated Fragment Width'}{'value'} = $cfragmentwidth; $Control_OvQual{'Estimated Fragment Width'}{'score'} = 2;
  $Control_OvQual{'Estimated Tag Length'}{'value'} = $cpredictedtaglength; $Control_OvQual{'Estimated Tag Length'}{'score'} = 2;
} #end if peaksxls
 
#Stitched regions & Superenhancers
if ($rose_stitched && $rose_superstitched){
  print "Processing ROSE counts ...";
  my $enhancers = `wc -l $rose_stitched | awk '{print \$1-6}'`; chop $enhancers;
  my $superenhancers = `wc -l $rose_superstitched | awk '{print \$1-6}'`; chop $superenhancers;
  unless ($enhancers > 0 ) { $enhancers = 0; }  #making sure enhancers  is numeric
  unless ($superenhancers > 0 ) { $superenhancers = 0; } #making sure superenhancers is numeric
  print OUT "Total number of Linear Stitched Peaks (Enhancers),", $enhancers,"\n";
  print OUT "Total number of SE-like Enriched Regions (Superenhancers),", $superenhancers,"\n";
  
  #QCdash
  $OVAL{$output_counter++} = 'Linear Stitched Peaks';
  $OVAL{$output_counter++} = 'SE-like Enriched Regions';

  $OvQual{'Linear Stitched Peaks'}{'value'} = $enhancers; $OvQual{'Linear Stitched Peaks'}{'score'} = -2;
  $OvQual{'SE-like Enriched Regions'}{'value'} = $superenhancers; $OvQual{'SE-like Enriched Regions'}{'score'} = 2;
  if ($OvQual{'Linear Stitched Peaks'}{'value'} > 1000) {$OvQual{'Linear Stitched Peaks'}{'score'} = -1;}
  if ($OvQual{'Linear Stitched Peaks'}{'value'} > 2000) {$OvQual{'Linear Stitched Peaks'}{'score'} = 0;}
  if ($OvQual{'Linear Stitched Peaks'}{'value'} > 5000) {$OvQual{'Linear Stitched Peaks'}{'score'} = 1;}
  if ($OvQual{'Linear Stitched Peaks'}{'value'} >= 10000) {$OvQual{'Linear Stitched Peaks'}{'score'} = 2;}

  if (($OvQual{'Linear Stitched Peaks'}{'value'} == 0) || ($OvQual{'SE-like Enriched Regions'}{'value'} == 0)) { $OvQual{'SE-like Enriched Regions'}{'score'} = -2 }
  else {
    if (($OvQual{'SE-like Enriched Regions'}{'value'}/$OvQual{'Linear Stitched Peaks'}{'value'}) > 0.02) {$OvQual{'SE-like Enriched Regions'}{'score'} = 1;}
    if (($OvQual{'SE-like Enriched Regions'}{'value'}/$OvQual{'Linear Stitched Peaks'}{'value'}) > 0.05) {$OvQual{'SE-like Enriched Regions'}{'score'} = 0;}
    if (($OvQual{'SE-like Enriched Regions'}{'value'}/$OvQual{'Linear Stitched Peaks'}{'value'}) > 0.1) {$OvQual{'SE-like Enriched Regions'}{'score'} = -1;}
    if (($OvQual{'SE-like Enriched Regions'}{'value'}/$OvQual{'Linear Stitched Peaks'}{'value'}) >= 0.2) {$OvQual{'SE-like Enriched Regions'}{'score'} = -2;}
    if ($OvQual{'SE-like Enriched Regions'}{'value'} <= 1) {$OvQual{'SE-like Enriched Regions'}{'score'} = -2;}
  }
} # end if rose

close(OUT);

#Processing QC html dashboard
my ($count,$totalscore) = (0,0); 
foreach (keys %OvQual){
  $count++;
  $totalscore += $OvQual{$_}{'score'};
}
#color names
my %color = ( "-2" => "#FF0000", "-1" => "#FF8C00", "0" => "#FFFF00", "1" => "#ADFF2F", "2" => "#008000" ); #red #orangered #yellow #greenyellow #green

my $OverallQuality = "POOR"; #($totalscore/$count);
my $color = $color{-2}; my $score = -2;
#if (length($OverallQuality) > 7) { $OverallQuality = sprintf ("%.4f", $OverallQuality) }
#if ((sprintf ("%.1f", ($totalscore/$count))) >= -1) { $color = $color{-1}; }
#if ((sprintf ("%.1f", ($totalscore/$count))) >= 0) { $color = $color{0}; }
#if ((sprintf ("%.1f", ($totalscore/$count))) >= 1) { $color = $color{1}; }
#if ((sprintf ("%.1f", ($totalscore/$count))) >= 2) { $color = $color{2}; }
if ((sprintf ("%.1f", ($totalscore/$count))) >= -1) { $color = $color{-1}; $score = -1; $OverallQuality = "BELOW_AVERAGE"; }
if ((sprintf ("%.1f", ($totalscore/$count))) >= 0) { $color = $color{0}; $score = 0; $OverallQuality = "AVERAGE";}
if ((sprintf ("%.1f", ($totalscore/$count))) >= 1) { $color = $color{1}; $score = 1; $OverallQuality = "GOOD";}
if ((sprintf ("%.1f", ($totalscore/$count))) >= 2) { $color = $color{2}; $score = 2; $OverallQuality = "EXCELLENT";}

my $htmlheader = "<table class='results'><tr><th>";
my $textheader = "Sample_Name";
my $samplehtmlvalues = "<tr><td><center>".$samplename."</center></td>";
my $controlhtmlvalues = "<tr><td><center>".$samplename."</center></td>";
my $sampletextvalues = $samplename;
my $controltextvalues = $samplename;

if ($cfastqczip) {
  $htmlheader .= "FASTQ";
  $textheader = "FASTQ";
  $samplehtmlvalues = "<tr><td><center>SAMPLE</center></td>";
  $controlhtmlvalues = "<tr><td><center>CONTROL</center></td>";
  $sampletextvalues = "SAMPLE";
  $controltextvalues = "CONTROL";
} else {
  $htmlheader .= "Sample Name";
}

open (CONFIG, ">$configfile");
print CONFIG "$samplename\tvalue\tscore\n";

#Adding Overall Quality
#if ($totaloutput) {
  $htmlheader .= "</th><th>"."Overall Quality";
  $textheader .= "\tOverall_Quality";
  if ($cfastqczip) { $samplehtmlvalues .= "<td bgcolor='$color' rowSpan='2'><center>".$OverallQuality."</center></td>"; }
  else { $samplehtmlvalues .= "<td bgcolor='$color'><center>".$OverallQuality."</center></td>";}
  $sampletextvalues .= "\t$OverallQuality";
  $controltextvalues .= "\t$OverallQuality";
  print CONFIG "Overall_Quality\t$OverallQuality\t$score\n";
#}

if ($cfastqczip) {
  foreach (sort {$a <=> $b} keys %CONTROL_OVAL){
    my $convertheader = $CONTROL_OVAL{$_};
    $convertheader =~ s/ /_/g; #change space to underscore for txt file
    $controlhtmlvalues .="<td bgcolor='".$color{$Control_OvQual{$CONTROL_OVAL{$_}}{'score'}}."'><center>".$Control_OvQual{$CONTROL_OVAL{$_}}{'value'}."</center></td>";
    $controltextvalues .= "\t$Control_OvQual{$CONTROL_OVAL{$_}}{'value'}";
    print CONFIG "control-$convertheader\t$Control_OvQual{$CONTROL_OVAL{$_}}{'value'}\t$Control_OvQual{$CONTROL_OVAL{$_}}{'score'}\n";
  }
  $controlhtmlvalues .= "</tr>";
}

foreach (sort {$a <=> $b} keys %OVAL){
  $htmlheader .= "</th><th>".$OVAL{$_};
  my $convertheader = $OVAL{$_};
  $convertheader =~ s/ /_/g; #change space to underscore for txt file
  $textheader .= "\t$convertheader";
  unless (exists $CONTROL_OVAL{$_}) {
    $samplehtmlvalues .= "<td bgcolor='".$color{$OvQual{$OVAL{$_}}{'score'}}."'";
    if ($cfastqczip) { $samplehtmlvalues .= " rowSpan='2'"; }
    $samplehtmlvalues .= "><center>".$OvQual{$OVAL{$_}}{'value'}."</center></td>";
    $controltextvalues .= "\t$OvQual{$OVAL{$_}}{'value'}";
  } else {
    $samplehtmlvalues .="<td bgcolor='".$color{$OvQual{$OVAL{$_}}{'score'}}."'><center>".$OvQual{$OVAL{$_}}{'value'}."</center></td>";
  }
  $sampletextvalues .= "\t$OvQual{$OVAL{$_}}{'value'}";
  print CONFIG "$convertheader\t$OvQual{$OVAL{$_}}{'value'}\t$OvQual{$OVAL{$_}}{'score'}\n";
}
$htmlheader .= "</th></tr>"; $samplehtmlvalues .= "</tr>";

close(CONFIG);

open (OUT2, ">$htmlfile"); #creating htmlfile
print OUT2 $htmlheader, "\n", $samplehtmlvalues;
if ($cfastqczip) { print OUT2 "\n", $controlhtmlvalues, "\n"; } 
close (OUT2);

open (OUT3, ">$textfile"); #creating htmlfile
print OUT3 $textheader, "\n", $sampletextvalues, "\n";
if ($cfastqczip) { print OUT3 $controltextvalues, "\n"; }
close (OUT3);
