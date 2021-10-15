#10_07_2021
Genome Reference FASTA
- Mouse mm9 reference genome from UCSC
https://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz

### Syntax
# wget https://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz;
# mkdir chromFA; cd chromFA; tar -xzvf ../chromFa.tar.gz;
# for chrom in $(ls -1 | grep -v "_random"); do cat $chrom >> mm9.fasta; done;


Genome Reference bowtie INDEX
Genome reference bowtie indexes were generated using Bowtie v1.2.3


Gene Annotation GTF
- downloaded from ucsc
# wget -O mm9_refGene.gtf.gz https://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/genes/mm9.refGene.gtf.gz


Genome Blacklist
- UHS blacklists from https://github.com/Boyle-Lab/Blacklist/blob/master/lists/Blacklist_v1/mm9-blacklist.bed.gz
# wget -O mm9.excludelist.bed.gz https://github.com/Boyle-Lab/Blacklist/raw/master/lists/Blacklist_v1/mm9-blacklist.bed.gz
