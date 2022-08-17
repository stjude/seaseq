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
    git checkout 3.0
    cd dnanexus
    reorg_id=$(dx build -f /seaseq-reorg | jq -r '.id')
    echo "Reorg applet ID: ${reorg_id}"
    sed -ibak "s/applet-Fx40j6091FfQp1P99p6b5k2x/${reorg_id}/" extras.json
    cd ..
    timestamp=$(date +%s)
    wget -nv https://github.com/dnanexus/dxCompiler/releases/download/2.9.0/dxCompiler-2.9.0.jar
    echo '434b515609123f1092453eac87984027  dxCompiler-2.9.0.jar' > dxCompiler-2.9.0.jar.md5
    md5sum -c dxCompiler-2.9.0.jar.md5

    if grep -F "control_sraid" /home/dnanexus/job_input.json; then
        SEASEQ="workflows/workflows/evaluatecontrolsrr.wdl"
        sed -i "s/cloud=false/cloud=true/" peaseq-control.wdl
        sed -i "s/cloud=false/cloud=true/" seaseq-control.wdl
    elif grep -F "sample_sraid" /home/dnanexus/job_input.json; then
        SEASEQ="workflows/workflows/evaluatesamplesrr.wdl"
        sed -i "s/cloud=false/cloud=true/" seaseq-case.wdl
        sed -i "s/cloud=false/cloud=true/" peaseq-case.wdl
    else 
        if grep -F "control_R2_fastq" /home/dnanexus/job_input.json; then
            SEASEQ="peaseq-control.wdl"
        elif grep -F "control_R1_fastq" /home/dnanexus/job_input.json; then
            SEASEQ="seaseq-control.wdl"
        elif grep -F "sample_R2_fastq" /home/dnanexus/job_input.json; then
            SEASEQ="peaseq-case.wdl"
        elif grep -F "sample_R1_fastq" /home/dnanexus/job_input.json; then
            SEASEQ="seaseq-case.wdl"
        fi
    fi
    
    echo $SEASEQ
    sed -i "s/cloud=false/cloud=true/" $SEASEQ
    
    sed -i "s/import \"..\/tasks\/util\.wdl/import \"\/home\/dnanexus\/seaseq\/workflows\/tasks\/util\.wdl/" workflows/workflows/visualization.wdl
    sed -i "s/import \"..\/tasks\/bedtools\.wdl/import \"\/home\/dnanexus\/seaseq\/workflows\/tasks\/bedtools\.wdl/" workflows/workflows/motifs.wdl
    sed -i "s/import \"..\/tasks\//import \"\/home\/dnanexus\/seaseq\/workflows\/tasks\//" workflows/workflows/mapping.wdl
    sed -i "s/import \"..\/../\//import \"\/home\/dnanexus\/seaseq\//" workflows/workflows/evaluatesamplesrr.wdl
    sed -i "s/import \"..\/../\//import \"\/home\/dnanexus\/seaseq\//" workflows/workflows/evaluatecontrolsrr.wdl

    # compile SEAseq to dxCompiler-2.9.0
    dx mkdir -p "${DX_PROJECT_CONTEXT_ID}":/app-$timestamp/
    wf_id=$(java -jar dxCompiler-2.9.0.jar compile $SEASEQ -project "${DX_PROJECT_CONTEXT_ID}" -folder /app-$timestamp -force -extras dnanexus/extras.json)
    echo "Workflow ID: ${wf_id}"

    # specify output folder and input json
    out_folder=$(dx describe --json "$DX_JOB_ID" | jq -r '.folder')

    echo "LAUNCHING WORKFLOW"
    jq 'walk( if type == "object" then . | del(.genome) else . end)' /home/dnanexus/updated_input.json | jq 'walk( if type == "object" then with_entries( .key |= if startswith("$dnanexus") or startswith("project") or startswith("id") then . else sub( "^"; "stage-common.") end)  else . end )' > /home/dnanexus/wf_input.json

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
