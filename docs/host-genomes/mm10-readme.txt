#10_07_2021
Genome Reference FASTA
- Mouse mm10 reference genome from UCSC
https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz

### Syntax
# wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz;
# mkdir chromFA; cd chromFA; tar -xzvf ../chromFa.tar.gz;
# for chrom in $(ls -1 | grep -v "_random" | grep -v "chrUn"); do cat $chrom >> mm10.fasta; done


Genome Reference bowtie INDEX
Genome reference bowtie indexes were generated using Bowtie v1.2.3


Gene Annotation GTF
- downloaded from ucsc
# wget -O mm10_refGene.gtf.gz https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.refGene.gtf.gz


Genome Blacklist
- UHS blacklists from https://github.com/Boyle-Lab/Blacklist/blob/master/lists/mm10-blacklist.v2.bed.gz
# wget -O mm10.excludelist.bed.gz https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz
