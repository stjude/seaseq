"""
    Applet to be called at the end of the SEASEQ workflow.
    Applet reorganizes all files from the defaulted output
        folder to prefered organization structure.
"""

import dxpy


@dxpy.entry_point("main")
def main(reorg_conf___=None, reorg_status___=None):  # pylint: disable=unused-argument

    # find the output stage of the current analysis
    analysis_id = dxpy.describe(dxpy.JOB_ID)["analysis"]
    stages = dxpy.describe(analysis_id)["stages"]

    # retrieve the dictionary containing outputs
    output_map = [
        x["execution"]["output"] for x in stages if x["id"] == "stage-outputs"
    ][0]
    folder_location = dxpy.describe(analysis_id)["folder"]

    project_container = dxpy.DXProject(dxpy.PROJECT_CONTEXT_ID)

    # move required outputfiles to their preferred permanent folders
    for file_identifiers in output_map.values():
        if isinstance(file_identifiers, (list, tuple)):
            for indvfile in file_identifiers:
                default_location = dxpy.describe(
                    dxpy.describe(indvfile["$dnanexus_link"])["createdBy"]["job"]
                )["runInput"]["default_location"]
                folder = folder_location + "/" + default_location
                project_container.new_folder(folder, parents=True)
        
                file_container = dxpy.describe(indvfile["$dnanexus_link"])["project"]
                file_object = dxpy.bindings.DXFile(
                    indvfile["$dnanexus_link"], project=file_container
                )
                if file_container == dxpy.PROJECT_CONTEXT_ID:
                    file_object.move(folder)
                else:
                    cloned_file = file_object.clone(  # pylint: disable=unused-variable
                        dxpy.PROJECT_CONTEXT_ID, folder=folder
                    )
        elif isinstance(file_identifiers, dict):
            if '$dnanexus_link' in file_identifiers:
                default_location = dxpy.describe(
                dxpy.describe(file_identifiers["$dnanexus_link"])["createdBy"]["job"])["runInput"]["default_location"]
                folder = folder_location + "/" + default_location
                project_container.new_folder(folder, parents=True)
        
                file_container = dxpy.describe(file_identifiers["$dnanexus_link"])["project"]
                file_object = dxpy.bindings.DXFile(
                    file_identifiers["$dnanexus_link"], project=file_container
                )
                if file_container == dxpy.PROJECT_CONTEXT_ID:
                    file_object.move(folder)
                else:
                    cloned_file = file_object.clone(
                        dxpy.PROJECT_CONTEXT_ID, folder=folder
                    )
