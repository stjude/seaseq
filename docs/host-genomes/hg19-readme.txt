#10_07_2021
Genome Reference FASTA
- Human hg19 reference genome from GENCODE Release 19  (GRCh37.p13): http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz. We excluded scaffolds, assembly patches and haplotypes.

### Syntax:
# wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz
# gunzip  GRCh37.p13.genome.fa.gz 
## selecting only reference chromosomes
# head -n 50 GRCh37.p13.genome.fa > hg19_chroms.fa


Genome Reference bowtie INDEX
Genome reference bowtie indexes were generated using Bowtie v1.2.3


Gene Annotation
- Human hg19 GTF comprehesive gene annotation from GENCODE Release 19 (GRCh37.p13): 
http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz.


Genome Blacklist
- UHS blacklists from hg19-blacklist.v2.bed.gz: https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg19-blacklist.v2.bed.gz
# wget -O hg19.excludelist.bed.gz https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg19-blacklist.v2.bed.gz
