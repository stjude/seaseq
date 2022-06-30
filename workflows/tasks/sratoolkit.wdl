version 1.0

task fastqdump {
    input {
        String? sra_id
        Boolean cloud=false
        Int memory_gb = 5
        Int max_retries = 1
        Int ncpu = 20
    }
    command {
        # configuration for sra-toolkit
        mkdir -p $PWD/ncbi $HOME/.ncbi
        echo """
        /LIBS/GUID = \"randomuuid\"
        /config/default = \"false\"
        /repository/user/ad/public/apps/file/volumes/flatAd = \".\"
        /repository/user/ad/public/apps/refseq/volumes/refseqAd = \".\"
        /repository/user/ad/public/apps/sra/volumes/sraAd = \".\"
        /repository/user/ad/public/apps/sraPileup/volumes/ad = \".\"
        /repository/user/ad/public/apps/sraRealign/volumes/ad = \".\"
        /repository/user/ad/public/apps/wgs/volumes/wgsAd = \".\"
        /repository/user/ad/public/root = \".\"
        /repository/user/default-path = \"$PWD/ncbi\"
        """ > $HOME/.ncbi/user-settings.mkfg

        # make number of threads compatible for DNAnexus
        threads=~{ncpu}
        if [ "~{cloud}" == 'true' ]; then
            output=$(sra-stat --meta --quick ~{sra_id})
            total=0
            for line in $output; do
                value=$(echo $line | cut -d '|' -f 3 | cut -d ':' -f 1)
                total=$(($total + $value))
            done
            if [ "$total" -ge 50000000 ]; then
                if [ "$threads" -ge 10 ]; then threads=10; fi
                echo "number of threads changed to 10 for dnanexus compatibility"
            fi
        fi

        # check if paired ended
        numLines=$(fastq-dump -X 1 -Z --split-spot ~{sra_id} | wc -l)
        paired_end="false"
        if [ "$numLines" -eq 8 ]; then
            paired_end="true"
        fi

        # perform fastqdump
        if [ "$paired_end" == 'true' ]; then
            echo true > paired_file
            pfastq-dump \
                -t $threads \
                --gzip \
                --split-files \
                -s ~{sra_id} -O ./

            #zcat -f ~{sra_id}_1.fastq.gz ~{sra_id}_2.fastq.gz | gzip -nc > ~{sra_id}.merged.fastq.gz
        else
            touch paired_file
            pfastq-dump \
                -t $threads \
                --gzip \
                -s ~{sra_id} -O ./
       fi 
    }
    runtime {
        memory: ceil(memory_gb * ncpu) + " GB"
        maxRetries: max_retries
        docker: 'ghcr.io/stjude/abralab/sratoolkit:v3.0.0'
        cpu: ncpu
    }
    output {
        File? R1end = "~{sra_id}_1.fastq.gz"
        File? R2end = "~{sra_id}_2.fastq.gz"
        Array[File] fastqfile = glob("~{sra_id}*fastq.gz")
        String paired_end = read_string('paired_file')
    }
}
