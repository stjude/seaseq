## STEPs to build DNAnexus workflow

1. Activate dx

    ```conda activate dx```

1. Login to account

    ```dx login```

1. Create New Project (if applicable)
    ```dx new project```

1. Build reorganization applet

    ```dx build seaseq-reorg```

1. Edit extras.json file with new applet id for seaseq-reorg

1. Build wdl workflow

    ```
    #using dxWDL-v1.50
    java -jar dxWDL-v1.50.jar compile seaseq.wdl \
        -project project-id \
        -extras extras.json \
        -folder /apps
    ```

1. Upload test data

    ``` dx upload -r test ```

1. Run workflow in interactive mode or with inputs specified

    ``` 
    #interactive mode
    dx run /apps/seaseq

    #or with the inputs specified
    dx run /apps/seaseq \
    -istage-common.blacklist="/test/hg19/hg19-blacklist.v2_chr21.bed" \
    -istage-common.fastqfile="/test/fastqfiles/H3K27Ac-AB5_R1_001.fastq.gz" \
    -istage-common.gtf="/test/hg19/GRCh37_latest_genomic_chr21.gtf" \
    -istage-common.reference="/test/hg19/CHR21INDEX/hg19_chr21.fa" \
    -istage-common.bowtie_index="/test/hg19/CHR21INDEX/hg19_chr21.1.ebwt" \
    -istage-common.bowtie_index="/test/hg19/CHR21INDEX/hg19_chr21.2.ebwt" \
    -istage-common.bowtie_index="/test/hg19/CHR21INDEX/hg19_chr21.3.ebwt" \
    -istage-common.bowtie_index="/test/hg19/CHR21INDEX/hg19_chr21.4.ebwt" \
    -istage-common.bowtie_index="/test/hg19/CHR21INDEX/hg19_chr21.rev.1.ebwt" \
    -istage-common.bowtie_index="/test/hg19/CHR21INDEX/hg19_chr21.rev.2.ebwt" \
    -istage-common.motif_databases="/test/hg19/motifs/CIS-BP.Homo_sapiens.meme" \
    -istage-common.motif_databases="/test/hg19/motifs/JASPAR2018_CORE_vertebrates_redundant.meme" \
    --destination "output-test"
    ```
