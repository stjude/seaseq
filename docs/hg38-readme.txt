#10_07_2021
Genome Reference FASTA
- Human hg38 reference genome from GENCODE Release 38  (GRCh38.p13): http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.p13.genome.fa.gz. We excluded scaffolds, assembly patches and haplotypes.

### Syntax:
# wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.p13.genome.fa.gz
# gunzip GRCh38.p13.genome.fa.gz
## selecting only reference chromosomes
# head -n 51471479 GRCh38.p13.genome.fa > hg38_chroms.fa


Genome Reference bowtie INDEX
Genome reference bowtie indexes were generated using Bowtie v1.2.3


Gene Annotation
- Human hg38 GTF comprehesive gene annotation from GENCODE Release 38 (GRCh38.p13): http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz.


Genome Blacklist
- UHS blacklists from hg38-blacklist.v2.bed.gz: https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
# wget -O hg38.excludelist.bed.gz https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
