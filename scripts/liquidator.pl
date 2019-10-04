#!/usr/bin/perl
# processing bamliquidator
use Pod::Usage;
use strict; 
use warnings;
use File::Basename;
use Getopt::Long;
use Statistics::R;

my $PATH = "/rgs01/project_space/abrahgrp/Software_Dev_Sandbox/common/ChipSeqPipeline/scripts";
my ($help, $manual, $rmdupbam, $gfffile, $outfile);
my $usage = "perl $0 -g <gff file> -b <bam file> [-outfile <outputfile>]\n";

GetOptions ("b|bam=s"=>\$rmdupbam,"g|gff=s"=>\$gfffile,"outfile|o=s"=>\$outfile);
unless ($rmdupbam && $gfffile) { die $usage; }
unless ($outfile) { $outfile = fileparse($rmdupbam, qr/\.[^.]*(\.bam)?$/); } 
else { $outfile = fileparse($outfile, qr/\.[^.]*(\..*)?$/); }

#generation of gff regions files.
`grep "gbkey=Gene" $gfffile | grep "gene_biotype=protein_coding" > genes.gff`;
`perl $PATH/flanking.pl -i genes.gff -f 2000 > promoters.gff`;
`bedtools flank -i genes.gff -g ~/.genomes/hg19/UCSC_CHROMSIZES/UCSC_hg19_chromInfo.tab -l 2000 -r 0 -s > upstream.gff`;
`bedtools flank -i genes.gff -g ~/.genomes/hg19/UCSC_CHROMSIZES/UCSC_hg19_chromInfo.tab -l 0 -r 2000 -s > downstream.gff`;

#creating bamliquidator files
print "running bam2GFF\n";
`python $PATH/bamToGFF/bamToGFF.py -b $rmdupbam -i ./promoters.gff -m 100 -o $outfile-promoters.txt`;
`python $PATH/bamToGFF/bamToGFF.py -b $rmdupbam -i ./upstream.gff -m 50 -o $outfile-upstream.txt`;
`python $PATH/bamToGFF/bamToGFF.py -b $rmdupbam -i ./downstream.gff -m 50 -o $outfile-downstream.txt`;
`python $PATH/bamToGFF/bamToGFF.py -b $rmdupbam -i ./genes.gff -m 100 -o $outfile-genebody.txt`;

#my @RCODE = qw|promoters.txt upstream.txt downstream.txt genebody.txt|;
#my @ROUTPUT = qw|promoters.pdf upstream.pdf downstream.pdf genebody.pdf|;
#my @RTITLE = qw|Bam-B100-Promoters Bam-B50-Upstream Bam-B50-Downstream Bam-B100-GeneBody|;
#my @COL = qw|red blue green black|;

#foreach (0..$#RCODE) { 
#  my $R = Statistics::R->new();
#  $R->startR;
#  $R->run("d <- read.table(\"$outfile-$RCODE[$_]\",sep=\"\t\",header=T);");
#  $R->run("pdf(\"$outfile-$ROUTPUT[$_]\")");
#  $R->run("matplot(colMeans(d[,3:ncol(d)]), type='l', col=\"$COL[$_]\", main=\"$RTITLE[$_]\")");
#  $R->run('dev.off()');
#  $R->stopR();
#}

my $R = Statistics::R->new();
$R->startR;
$R->run("promoters <- read.table(\"$outfile-promoters.txt\",sep=\"\\t\",header=T);");
$R->run("upstream <- read.table(\"$outfile-upstream.txt\",sep=\"\\t\",header=T);");
$R->run("downstream <- read.table(\"$outfile-downstream.txt\",sep=\"\\t\",header=T);");
$R->run("genebody <- read.table(\"$outfile-genebody.txt\",sep=\"\\t\",header=T);");
$R->run("pdf(\"$outfile-promoters.pdf\")");
$R->run("matplot(colMeans(promoters[,3:ncol(promoters)]), type='l', main=\"Bam-B100-Promoters\")");
$R->run('dev.off()');
$R->run("pdf(\"$outfile-entiregene.pdf\")");
$R->run('combined<-cbind(upstream[,3:ncol(upstream)], genebody[,3:ncol(genebody)], downstream[,3:ncol(downstream)])');
$R->run('matplot(colMeans(combined),type="l",main="Bam-MetaGene")');
$R->run('dev.off()');
$R->stopR();
