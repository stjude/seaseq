#!/bin/bash

set -e -o -x

main() {
    git clone https://github.com/stjude/seaseq.git seaseq
    cd seaseq
    git checkout 1.0
    cd dnanexus
    reorg_id=$(dx build -f /seaseq-reorg | jq -r '.id')
    echo "Reorg applet ID: ${reorg_id}"
    sed -ibak "s/applet-Fx40j6091FfQp1P99p6b5k2x/${reorg_id}/" extras.json
    cd ..

    wget -nv https://github.com/dnanexus/dxWDL/releases/download/v1.50/dxWDL-v1.50.jar
    echo '476621564b3b310b17598ee1f02a1865 dxWDL-v1.50.jar' > dxWDL-v1.50.jar.md5
    md5sum -c dxWDL-v1.50.jar.md5

    sed -ibak "s/cloud=\"false\"/cloud=\"true\"/" seaseq.wdl
    dx mkdir -p "${DX_PROJECT_CONTEXT_ID}":/apps/
    wf_id=$(java -jar dxWDL-v1.50.jar compile seaseq.wdl -project "${DX_PROJECT_CONTEXT_ID}" -folder /apps -force -extras dnanexus/extras.json)
    echo "Workflow ID: ${wf_id}"

    out_folder=$(dx describe --json "$DX_JOB_ID" | jq -r '.folder')
    echo "LAUNCHING WORKFLOW"
    jq 'walk( if type == "object" then with_entries( .key |= if startswith("$dnanexus") then . else sub( "^"; "stage-common.") end)  else . end )' /home/dnanexus/job_input.json > /home/dnanexus/wf_input.json
    analysis_id=$(dx run -y "$wf_id" --brief -f /home/dnanexus/wf_input.json --destination "${DX_PROJECT_CONTEXT_ID}":"$out_folder" --extra-args '{"executionPolicy": {"onNonRestartableFailure": "failStage"}}')
    dx wait "$analysis_id"
    dx describe --json "$analysis_id" \
        | jq ".output" \
        | grep -v "reorg_status" \
        | jq 'walk(if type == "object"  then with_entries(select(.key | test("^stage-outputs|dnanexus_link|^___$") )) else . end) | walk( if type == "object" then with_entries( .key |= if startswith("stage-outputs") then sub( "^stage-outputs."; "") else . end  )  else . end ) |walk( if type == "object" then with_entries( .key |= if endswith("___dxfiles") then sub( "___dxfiles$"; "") else . end  )  else . end )  | walk (if (type == "object" and has("___")) then .|=.___ else . end)' > /home/dnanexus/job_output.json

    dx rm -rf "${DX_PROJECT_CONTEXT_ID}":/apps/
}
