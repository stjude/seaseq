#!/usr/bin/perl
#normalize peak calls

$inputWIG = $ARGV[0]; #input wig file
$xlsfile = $ARGV[1]; #peak calls summary file
$inputWIG =~ /(.*).wig/;
$outputWIG = ((split("/",$1))[-1]).'.RPM.wig'; #generate output file
$mappedReads = `grep "tags after filtering in treatment" $xlsfile | awk -F':' '{print \$2}'`;
$mappedReads =~ s/\s//g;
unless ($mappedReads > 0) { #taking total tags in treatment instead
  $mappedReads = `grep "total tags in treatment" $xlsfile | awk -F':' '{print \$2}'`;
  $mappedReads =~ s/\s//g;
} 

$mappedReads = $mappedReads/1000000; #per million mapped reads.

open(FD,"gunzip -c $inputWIG |");
open(FOUT,"> $outputWIG");
print "done";
while(<FD>){
  chomp;
  if($_ !~ m/track/ && $_ !~ m/variable/) {
    @lineArr = split("\t", $_);
    my $height = $lineArr[1];
    $height = $height/$mappedReads;
    print FOUT $lineArr[0] . "\t" . $height."\n";
  } else{
    my $header = ($_ =~ s/name\=\"/name\=\"RPM\_/);
    print FOUT $_."\n";
  }
}
close(FD);
close(FOUT);
`gzip $outputWIG`;
