## STEPs to build DNAnexus workflow
```
conda activate dx

#login to account
dx login

#create new project if applicable
dx new project

#build reorganization applet
dx build seaseq-reorg

#build wdl workflow
java -jar dxWDL-v1.48.2.jar compile seaseq.wdl \
    -project project-id \
    -reorg -extras.json \
    -folder /apps

#upload test data
dx upload -r test

#run workflow in interactive mode
dx run /apps/seaseq

#or with the inputs specified
dx run /apps/seaseq \
-istage-common.blacklistfile="/test/hg19/hg19-blacklist.v2_chr21.bed" \
-istage-common.chromsizes="/test/hg19/UCSC_hg19_chromInfo_chr21.tab" \
-istage-common.fastqfile="/test/fastqfiles/H3K27Ac-AB5_R1_001.fastq.gz" \
-istage-common.gtffile="/test/hg19/GRCh37_latest_genomic_chr21.gtf" \
-istage-common.reference="/test/hg19/CHR21INDEX/hg19_chr21.fa" \
-istage-common.reference_index="/test/hg19/CHR21INDEX/hg19_chr21.fa.fai" \
-istage-common.index_files="/test/hg19/CHR21INDEX/hg19_chr21.1.ebwt" \
-istage-common.index_files="/test/hg19/CHR21INDEX/hg19_chr21.2.ebwt" \
-istage-common.index_files="/test/hg19/CHR21INDEX/hg19_chr21.3.ebwt" \
-istage-common.index_files="/test/hg19/CHR21INDEX/hg19_chr21.4.ebwt" \
-istage-common.index_files="/test/hg19/CHR21INDEX/hg19_chr21.rev.1.ebwt" \
-istage-common.index_files="/test/hg19/CHR21INDEX/hg19_chr21.rev.2.ebwt" \
-istage-common.motif_databases="/test/hg19/motifs/CIS-BP.Homo_sapiens.meme" \
-istage-common.motif_databases="/test/hg19/motifs/JASPAR2018_CORE_vertebrates_redundant.meme" \
--destination "output-test"

```
