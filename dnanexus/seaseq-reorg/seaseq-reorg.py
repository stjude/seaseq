'''
    Applet to be called at the end of the SEASEQ workflow.
    Applet reorganizes all files from the defaulted output
        folder to prefered organization structure.
'''

from datetime import datetime

import dxpy

@dxpy.entry_point('main')

def main(reorg_conf___=None, reorg_status___=None):

    # find the output stage of the current analysis
    analysis_id = dxpy.describe(dxpy.JOB_ID)["analysis"]
    stages = dxpy.describe(analysis_id)["stages"]

    # retrieve the dictionary containing outputs
    output_map = [x['execution']['output'] for x in stages if x['id'] == 'stage-outputs'][0]
    common_map =  [x['execution']['output'] for x in stages if x['id'] == 'stage-common'][0]

    folder_location = [x['execution']['folder'] for x in stages if x['id'] == 'stage-outputs'][0]

    # retrieve container id
    dx_container = dxpy.DXProject(dxpy.PROJECT_CONTEXT_ID)

    # create temporary folder to move all created files in the output folder
    datestamp = datetime.today().strftime('%s')
    temp_folder = "/" + datestamp + "-temp"
    dx_container.new_folder(temp_folder, parents=True)

    #pipeline input files are moved to the temporary folder if applicable
    for eachcommon in common_map.values():
        if isinstance(eachcommon, dict):
            common_folder = dxpy.describe(eachcommon['$dnanexus_link'])['folder']
            if folder_location == common_folder:
                dx_container.move(
                    destination=temp_folder,
                    objects=[eachcommon['$dnanexus_link']]
                )
        elif isinstance(eachcommon,list):
            for thepath in eachcommon:
                common_folder = dxpy.describe(thepath['$dnanexus_link'])
                if folder_location == common_folder:
                    dx_container.move(
                        destination=temp_folder,
                        objects=[thepath['$dnanexus_link']]
                    )

    #intermediate files. All files are moved to the temporary folder
    for eachfile in dx_container.list_folder(folder_location)['objects']:
        checkjobid = dxpy.describe(eachfile['id'])['createdBy']
        if 'job' in checkjobid.keys():
            if dxpy.describe(checkjobid['job'])['rootExecution'] == analysis_id:
                dx_container.move(
                    destination=temp_folder,
                    objects=[eachfile['id']]
                )

    # move required outputfiles to their preferred permanent folders
    for file_identifiers in output_map:
        file_id_information = output_map[file_identifiers]
        if isinstance(file_id_information, list):
            for indvfile in range(0,len(file_id_information)):
                default_location = dxpy.describe(dxpy.describe(file_id_information[indvfile])['createdBy']['job'])['runInput']['default_location']
                folder = folder_location + "/" + default_location
                dx_container.new_folder(folder, parents=True)
                dx_container.move(
                    destination=folder,
                    objects=[file_id_information[indvfile]['$dnanexus_link']]
                )
                
    # remove the temporary folder created
    dx_container.remove_folder(folder=temp_folder, recurse=True)
