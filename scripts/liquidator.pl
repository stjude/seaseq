#!/usr/bin/perl
# processing bamliquidator
use Pod::Usage;
use strict; 
use warnings;
use File::Basename;
use Getopt::Long;

my $PATH = "/rgs01/project_space/abrahgrp/Software_Dev_Sandbox/common/madetunj/software";
my ($help, $manual, $rmdupbam, $gtffile, $outfile, $samplename);
my (%HASH, %CONTENT);
my $usage = "perl $0 -g <gtf file> -b <bam file> [-outfile <outputfile>] [-sample <samplename>]\n";

GetOptions ("b|bam=s"=>\$rmdupbam,"g|gtf=s"=>\$gtffile,"outfile|o=s"=>\$outfile, "sample|s=s"=>\$samplename);
unless ($rmdupbam && $gtffile) { die $usage; }
unless ($outfile) { $outfile = fileparse($rmdupbam, qr/(\.bam)?$/); } 
else { $outfile = fileparse($outfile, qr/(\.[\d\w]+)?$/); }
unless ($samplename) {$samplename = fileparse($outfile, qr/\.[^.]*(\..*)?$/); }

#generation of gff file
#initially gff file 
#`grep "gbkey=Gene" $gfffile | grep "gene_biotype=protein_coding" > genes.gff`;
#now from gtf file
open (GTF, "<", $gtffile); open (GFF, '>genes.gff');
while (<GTF>) {
  unless (/^#/) {
    chomp;
    if ($_ =~ /"protein_coding"/) {
      $_ =~ s/\"//g;
      my @l = split /\t/; 
      if ($l[2] =~ /^gene$/) {
        $l[0] =~ s/^M.+/M/g; #changing MT to M
        print GFF "chr",$l[0],"\t"; #added chr to chrom
        print GFF join("\t",@l[1..7]); #plus the 2nd to 8th column
        my @line = split(" ", $l[8]); #split 8th column
        print GFF "\t",$line[0],"=",$line[1],$line[4],"=",$line[5],$line[8],"=",$line[9],"\n"; #plus geneid,genename,genebiotype
      }
    }
  }
}
close (GTF); close (GFF);

#generating the gff regions files.
`flanking.pl -i genes.gff -f 2000 > promoters.gff`;
`bedtools flank -i genes.gff -g ~/.genomes/hg19/UCSC_CHROMSIZES/UCSC_hg19_chromInfo.tab -l 2000 -r 0 -s > upstream.gff`;
`bedtools flank -i genes.gff -g ~/.genomes/hg19/UCSC_CHROMSIZES/UCSC_hg19_chromInfo.tab -l 0 -r 2000 -s > downstream.gff`;

#creating bamliquidator files
print "running bam2GFF\n";
`python $PATH/bamToGFF/bamToGFF.py -b $rmdupbam -i ./promoters.gff -m 100 -o $outfile-promoters.txt`;
`python $PATH/bamToGFF/bamToGFF.py -b $rmdupbam -i ./upstream.gff -m 50 -o $outfile-upstream.txt`;
`python $PATH/bamToGFF/bamToGFF.py -b $rmdupbam -i ./downstream.gff -m 50 -o $outfile-downstream.txt`;
`python $PATH/bamToGFF/bamToGFF.py -b $rmdupbam -i ./genes.gff -m 100 -o $outfile-genebody.txt`;

foreach my $index (qw|promoters upstream downstream genebody|){
  open (IN, "<$outfile-$index.txt");
  my $linenumber = 0; 
  while (<IN>) {
    chomp;
    my @line = split /\t/;
    unless ($line[0] =~ /GENE_ID/) {
      $linenumber++;
      my $sum = 0; 
      foreach my $i (@line[2..$#line]){ if ($i =~ /NA/) { $i = 0; } $sum = $sum + $i; } 
      $HASH{$sum}{$linenumber} = $sum;
      $CONTENT{$linenumber} = $_; 
    } 
  }
  close (IN);
  open(OUT,">$outfile-sorted.$index.txt");
  foreach my $asum (sort {$b <=> $a} keys %HASH) {
    foreach my $aline (keys %{$HASH{$asum}}) {
      print OUT $CONTENT{$aline},"\n";
      delete $HASH{$asum}{$aline};
      delete $CONTENT{$aline};
    }
  }
  close (OUT);
}

my $Rcode = <<"ENDOFR";
#!/usr/bin/env Rscript
source(\"$PATH/heatmap.3.R\");
#read in original files
promoters <- read.table("$outfile-sorted.promoters.txt",sep="\\t",header=F);
upstream <- read.table("$outfile-sorted.upstream.txt",sep="\\t",header=F);
downstream <- read.table("$outfile-sorted.downstream.txt",sep="\\t",header=F);
genebody <- read.table("$outfile-sorted.genebody.txt",sep="\\t",header=F);

#combining entire genebody
combined<-cbind(upstream[,3:ncol(upstream)], genebody[,3:ncol(genebody)], downstream[,3:ncol(downstream)]);

#removing NA.
promoters<-na.omit(promoters)
combined<-na.omit(combined)
#matplot of promoters & genebody

pdf("$outfile-promoters.pdf");
matplot(colMeans(promoters[,3:ncol(promoters)]), type='l', main="$samplename Promoters", ylab="Average normalized mapped reads", xlim=NULL, xaxt='n', xlab="Genomic Region");
axis(1, at=c(0,50,100), labels=c("-50", "TSS", "50"));
dev.off();
pdf("$outfile-entiregene.pdf");
matplot(colMeans(combined),type='l',main="$samplename  MetaGenes", ylab="Average normalized mapped reads", xlim=NULL, xaxt='n', xlab="Genomic Region");
axis(1, at=c(0,50,83,116,150,200), labels=c("-50", "TSS", "33%","66%", "TES", "50"));
dev.off();

#heatmap of promoters & genebody
colz=colorRampPalette(c("white", "red"))(quantile(as.vector(t(promoters[,3:ncol(promoters)])),.80));
breaks=seq(0,(quantile(as.vector(t(promoters[,3:ncol(promoters)])),.80))+1,by=1);
png("$outfile-heatmap.promoters.png", type="cairo");
heatmap.3(promoters[,3:ncol(promoters)], col=colz, breaks=breaks, trace="none", dendrogram="none", Colv=NA, Rowv=NA, density.info="none", labRow=NA, labCol=NA, main="$samplename\nPromoters");
dev.off();
colz=colorRampPalette(c("white", "red"))(quantile(as.vector(t(combined[,3:ncol(combined)])),.80));
breaks=seq(0,(quantile(as.vector(t(combined[,3:ncol(combined)])),.80))+1,by=1);
png("$outfile-heatmap.entiregene.png",type="cairo");
heatmap.3(combined[,3:ncol(combined)], col=colz, breaks=breaks, trace="none", dendrogram="none", Colv=NA, Rowv=NA, density.info="none", labRow=NA, labCol=NA, main="$samplename\nMetaGenes");
dev.off();
ENDOFR

open (OUTR, ">Rcodeout.R");
print OUTR $Rcode."\n";
close OUTR;

`Rscript ./Rcodeout.R`;
