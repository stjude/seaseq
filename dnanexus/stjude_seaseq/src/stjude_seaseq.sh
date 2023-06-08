#!/bin/bash

set -e -o -x

main() {
    # Check genome provided, either choosing from the genomes provided from the dropdown (genome)
    # or specify the genome files after choosing Custom
    # Make sure if a host genome or custom with required genome files are provided.
    echo "Value of genome: '$genome'"
    if [ "$genome" != "" -a ${genome:0:6} != "Custom" ]
    then
      if [ "$reference" != "" ]
      then
        echo "Could not determine which Genome to use." >&2
        echo "A custom Genome Reference was provided but Host Genome was not set to custom." >&2
        exit 1
      fi
      if [ "$gtf" != "" ]
      then
        echo "Could not determine which Genome to use." >&2
        echo "A custom Gene Annotation was provided but Host Genome was not set to custom." >&2
        exit 1
      fi
      if [ "$bowtie_index" != "" ]
      then
        echo "Could not determine which Genome to use." >&2
        echo "A custom Bowtie Index was provided but Host Genome was not set to custom." >&2
        exit 1
      fi
      #adding genome config to the input.json
      jq -s add /genomes/${genome##*/}-config.json /home/dnanexus/job_input.json > /home/dnanexus/updated_input.json

    else
      if [ "$reference" == "" ]
      then
        echo "Could not determine which Genome to use." >&2
        echo "A custom Genome Reference was not provided but Host Genome was set to custom." >&2
        exit 1
      fi
      if [ "$gtf" == "" ]
      then
        echo "Could not determine which Genome to use." >&2
        echo "A custom Gene Annotation was not provided but Host Genome was set to custom." >&2
        exit 1
      fi
      cp /home/dnanexus/job_input.json /home/dnanexus/updated_input.json

    fi

    # Make sure sample inputs are provided
    if [ "$sample_R1_fastq" == "" -a "$sample_sraid" == "" ]
    then
      echo "Sample Input could not be determined." >&2
      echo "Sample FASTQ or SRR was not provided." >&2
      exit 1
    fi

    # Make sure required files are provided
    if grep -F "reference" /home/dnanexus/updated_input.json; then echo "STATUS : Genome FASTA provided"
    else
      echo "Genome Reference FASTA not provided" >&2
      exit 1
    fi
    if grep -F "gtf" /home/dnanexus/updated_input.json; then echo "STATUS : Gene Annotation provided"
    else
      echo "Gene Annotation not provided" >&2
      exit 1
    fi

    # git clone and configure SEAseq
    git clone https://github.com/stjude/seaseq.git seaseq
    cd seaseq
    git checkout 3.1
    cd dnanexus
    reorg_id=$(dx build -f /seaseq-reorg | jq -r '.id')
    echo "Reorg applet ID: ${reorg_id}"
    sed -ibak "s/applet-Fx40j6091FfQp1P99p6b5k2x/${reorg_id}/" extras.json
    cd ..
    timestamp=$(date +%s)
    #wget -nv https://github.com/dnanexus/dxWDL/releases/download/v1.50/dxWDL-v1.50.jar
    #echo '476621564b3b310b17598ee1f02a1865 dxWDL-v1.50.jar' > dxWDL-v1.50.jar.md5
    #md5sum -c dxWDL-v1.50.jar.md5

    wget -nv https://github.com/dnanexus/dxCompiler/releases/download/2.11.2/dxCompiler-2.11.2.jar
    echo '64c05d41bed875317b7e6cfe1db588b1 dxCompiler-2.11.2.jar' > dxCompiler-2.11.2.jar.md5
    md5sum -c dxCompiler-2.11.2.jar.md5

    # Check fastqtype
    if [ $fastqtype = true ]; then pipeline_prefix='peaseq';
    else
        pipeline_prefix='seaseq';
    fi

    # Final pipeline selection
    if grep -F "control_sraid" /home/dnanexus/job_input.json; then
        SEASEQ=$pipeline_prefix"-control.wdl"
    elif grep -F "sample_sraid" /home/dnanexus/job_input.json; then
        if grep -F "control_R1_fastq" /home/dnanexus/job_input.json; then
            SEASEQ=$pipeline_prefix"-control.wdl"
        else
            SEASEQ=$pipeline_prefix"-case.wdl"
        fi
    else 
        if grep -F "control_R2_fastq" /home/dnanexus/job_input.json; then
            SEASEQ="peaseq-control.wdl"
        elif grep -F "control_R1_fastq" /home/dnanexus/job_input.json; then
            SEASEQ="seaseq-control.wdl"
            sed -i "s/sample_R1_fastq/sample_fastq/" /home/dnanexus/updated_input.json
            sed -i "s/control_R1_fastq/control_fastq/" /home/dnanexus/updated_input.json
        elif grep -F "sample_R2_fastq" /home/dnanexus/job_input.json; then
            SEASEQ="peaseq-case.wdl"
        elif grep -F "sample_R1_fastq" /home/dnanexus/updated_input.json; then
            SEASEQ="seaseq-case.wdl"
            sed -i "s/sample_R1_fastq/sample_fastq/" /home/dnanexus/updated_input.json
        else
            echo "Can Not determine Pipeline to execute, contact Maintainer" >&2
            exit 1
        fi
    fi

    echo $SEASEQ;
    if [[ $SEASEQ == "seaseq-case.wdl" || $SEASEQ == "seaseq-control.wdl" ]]; then
        jq 'walk( if type == "object" then . | del(.insertsize) else . end)' /home/dnanexus/updated_input.json | jq 'walk( if type == "object" then . | del(.strandedness) 
else . end)' > /home/dnanexus/corrected_input.json
        echo "yes"
    else
        cp /home/dnanexus/updated_input.json /home/dnanexus/corrected_input.json
        echo "no"
    fi

    sed -i "s/cloud=false/cloud=true/" $SEASEQ
    grep cloud $SEASEQ
    before='import "../tasks/'; before="${before//\//\\/}"
    after='import "/home/dnanexus/seaseq/workflows/tasks/'; after="${after//\//\\/}"
    sed -i "s/${before}/${after}/" workflows/workflows/visualization.wdl
    sed -i "s/${before}/${after}/" workflows/workflows/motifs.wdl
    sed -i "s/${before}/${after}/" workflows/workflows/mapping.wdl
    
    # compile SEAseq to dxWDL-v1.50 
    dx mkdir -p "${DX_PROJECT_CONTEXT_ID}":/app-$timestamp/
    #wf_id=$(java -jar dxWDL-v1.50.jar compile $SEASEQ -project "${DX_PROJECT_CONTEXT_ID}" -folder /app-$timestamp -force -extras dnanexus/extras.json)
    
    wf_id=$(java -jar dxCompiler-2.11.2.jar compile $SEASEQ -project "${DX_PROJECT_CONTEXT_ID}" -folder /app-$timestamp -force -extras dnanexus/extras.json)
    echo "Workflow ID: ${wf_id}"

    # specify output folder and input json
    out_folder=$(dx describe --json "$DX_JOB_ID" | jq -r '.folder')

    echo "LAUNCHING WORKFLOW"
    jq 'walk( if type == "object" then . | del(.genome) else . end)' /home/dnanexus/corrected_input.json | jq 'walk( if type == "object" then . | del(.fastqtype) else . end)' | jq 'walk( if type == "object" then with_entries( .key |= if startswith("$dnanexus") or startswith("project") or startswith("id") then . else sub( "^"; "stage-common.") end)  else . end )' > /home/dnanexus/wf_input.json

    cat /home/dnanexus/wf_input.json

    analysis_id=$(dx run -y "$wf_id" --brief -f /home/dnanexus/wf_input.json --destination "${DX_PROJECT_CONTEXT_ID}":"$out_folder" --extra-args '{"executionPolicy": {"onNonRestartableFailure": "failStage"}}')

    # LAUNCH SEAseq
    dx wait "$analysis_id"

    # wait until analysis is completed and export output files
    dx describe --json "$analysis_id" \
        | jq ".output" \
        | grep -v "reorg_status" \
        | jq 'walk(if type == "object"  then with_entries(select(.key | test("^stage-outputs|dnanexus_link|^___$") )) else . end) | walk( if type == "object" then with_entries( .key |= if startswith("stage-outputs") then sub( "^stage-outputs."; "") else . end  )  else . end ) |walk( if type == "object" then with_entries( .key |= if endswith("___dxfiles") then sub( "___dxfiles$"; "") else . end  )  else . end )  | walk (if (type == "object" and has("___")) then .|=.___ else . end)' > /home/dnanexus/job_output.json

    # remove configured SEAseq
    dx rm -rf "${DX_PROJECT_CONTEXT_ID}":/app-$timestamp/
}

