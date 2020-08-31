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
    folder_location = [x['execution']['folder'] for x in stages if x['id'] == 'stage-outputs'][0]

    # retrieve container id
    dx_container = dxpy.DXProject(dxpy.PROJECT_CONTEXT_ID)

    # create temporary folder to move all files in the output folder
    datestamp = datetime.today().strftime('%s')
    temp_folder = "/" + datestamp + "-temp"
    dx_container.new_folder(temp_folder, parents=True)
    for eachfile in dx_container.list_folder(folder_location)['objects']:
        dx_container.move(
            destination=temp_folder,
            objects=[eachfile['id']]
        )

    # move required outputfiles to their preferred folders
#FASTQC
    htmlfile = output_map['htmlfile']
    zipfile = output_map['zipfile']
    bam_htmlfile = output_map['bam_htmlfile']
    bam_zipfile = output_map['bam_zipfile']

    sample_name = dxpy.describe(htmlfile)["name"].split('_fastqc.html')[0]

    folder = folder_location + "/" + sample_name + "/QC_files/FASTQC"
    dx_container.new_folder(folder, parents=True)
    for thestring in (htmlfile, zipfile, bam_htmlfile, bam_zipfile):
        if thestring is not None:
            dx_container.move(
                destination=folder,
                objects=[thestring['$dnanexus_link']]
            )

#BASICMETRICS
    metrics_out = output_map['metrics_out']

#QC-STATS
    qc_statsfile = output_map['qc_statsfile']
    qc_htmlfile = output_map['qc_htmlfile']
    qc_textfile = output_map['qc_textfile']

    folder = folder_location + "/" + sample_name + "/QC_files/STATS"
    dx_container.new_folder(folder, parents=True)
    for thestring in (metrics_out, qc_statsfile, qc_htmlfile, qc_textfile):
        if thestring is not None:
            dx_container.move(
                destination=folder,
                objects=[thestring['$dnanexus_link']]
            )

#BAMFILES
    sortedbam = output_map['sortedbam']
    mkdupbam = output_map['mkdupbam']
    bklistbam = output_map['bklistbam']
    indexbam = output_map['indexbam']
    bklist_indexbam = output_map['bklist_indexbam']
    mkdup_indexbam = output_map['mkdup_indexbam']

    folder = folder_location + "/" + sample_name + "/BAM_files"
    dx_container.new_folder(folder, parents=True)
    for thestring in (sortedbam, mkdupbam, bklistbam, indexbam, bklist_indexbam, mkdup_indexbam):
        if thestring is not None:
            dx_container.move(
                destination=folder,
                objects=[thestring['$dnanexus_link']]
            )

#MACS
    peakbedfile = output_map['peakbedfile']
    peakxlsfile = output_map['peakxlsfile']
    summitsfile = output_map['summitsfile']
    wigfile = output_map['wigfile']
    all_peakbedfile = output_map['all_peakbedfile']
    all_peakxlsfile = output_map['all_peakxlsfile']
    all_summitsfile = output_map['all_summitsfile']
    all_wigfile = output_map['all_wigfile']
    nm_peakbedfile = output_map['nm_peakbedfile']
    nm_peakxlsfile = output_map['nm_peakxlsfile']
    nm_summitsfile = output_map['nm_summitsfile']
    nm_wigfile = output_map['nm_wigfile']

    prefix_peaks_folder = folder_location + "/" + sample_name + "/PEAKS_files"
    folder = prefix_peaks_folder + "/NARROW_peaks"
    dx_container.new_folder(folder, parents=True)
    for thestring in (peakbedfile, peakxlsfile, summitsfile, wigfile,
                      all_peakbedfile, all_peakxlsfile, all_summitsfile, all_wigfile,
                      nm_peakbedfile, nm_peakxlsfile, nm_summitsfile, nm_wigfile):
        if thestring is not None:
            dx_container.move(
                destination=folder,
                objects=[thestring['$dnanexus_link']]
            )

#SICER
    scoreisland = output_map['scoreisland']
    sicer_wigfile = output_map['sicer_wigfile']

    folder = prefix_peaks_folder + "/BROAD_peaks"
    dx_container.new_folder(folder, parents=True)
    for thestring in (scoreisland, sicer_wigfile):
        if thestring is not None:
            dx_container.move(
                destination=folder,
                objects=[thestring['$dnanexus_link']]
            )

#ROSE
    pngfile = output_map['pngfile']
    mapped_union = output_map['mapped_union']
    mapped_stitch = output_map['mapped_stitch']
    enhancers = output_map['enhancers']
    super_enhancers = output_map['super_enhancers']
    gff_file = output_map['gff_file']
    gff_union = output_map['gff_union']
    union_enhancers = output_map['union_enhancers']
    stitch_enhancers = output_map['stitch_enhancers']
    e_to_g_enhancers = output_map['e_to_g_enhancers']
    g_to_e_enhancers = output_map['g_to_e_enhancers']
    e_to_g_super_enhancers = output_map['e_to_g_super_enhancers']
    g_to_e_super_enhancers = output_map['g_to_e_super_enhancers']

    prefolder = prefix_peaks_folder + "/STITCHED_REGIONS"
    folder = prefolder + '/gff'
    dx_container.new_folder(folder, parents=True)
    for thestring in (gff_file, gff_union):
        if thestring is not None:
            dx_container.move(
                destination=folder,
                objects=[thestring['$dnanexus_link']]
            )
    folder = prefolder + '/mappedGFF'
    dx_container.new_folder(folder, parents=True)
    for thestring in (mapped_union, mapped_stitch):
        if thestring is not None:
            dx_container.move(
                destination=folder,
                objects=[thestring['$dnanexus_link']]
            )
    for thestring in (pngfile, enhancers, super_enhancers, union_enhancers,
                      stitch_enhancers, e_to_g_enhancers, e_to_g_super_enhancers,
                      g_to_e_enhancers, g_to_e_super_enhancers):
        if thestring is not None:
            dx_container.move(
                destination=prefolder,
                objects=[thestring['$dnanexus_link']]
            )

#FLANKBED
    flankbedfile = output_map['flankbedfile']

#MOTIFS
    ame_tsv = output_map['ame_tsv']
    ame_html = output_map['ame_html']
    ame_seq = output_map['ame_seq']
    meme = output_map['meme']
    meme_summary = output_map['meme_summary']
    summit_ame_tsv = output_map['summit_ame_tsv']
    summit_ame_html = output_map['summit_ame_html']
    summit_ame_seq = output_map['summit_ame_seq']
    summit_meme = output_map['summit_meme']
    summit_meme_summary = output_map['summit_meme_summary']

    prefolder = folder_location + "/" + sample_name + "/MOTIF_files"
    summit_motif_subfile = 'bklist' + \
                           (dxpy.describe(flankbedfile)["name"].split('bklist')[1]).split('.bed')[0]
    motif_subfile = 'bklist' + \
                    (dxpy.describe(peakbedfile)["name"].split('bklist')[1]).split('.bed')[0]
    folder = prefolder + "/" + motif_subfile + '-ame_out'
    dx_container.new_folder(folder, parents=True)
    for thestring in (ame_tsv, ame_html, ame_seq):
        if thestring is not None:
            dx_container.move(
                destination=folder,
                objects=[thestring['$dnanexus_link']]
            )
    folder = prefolder + "/" + summit_motif_subfile + '-ame_out'
    dx_container.new_folder(folder, parents=True)
    for thestring in (summit_ame_tsv, summit_ame_html, summit_ame_seq):
        if thestring is not None:
            dx_container.move(
                destination=folder,
                objects=[thestring['$dnanexus_link']]
            )
    #move flankbedfile and meme files to main MOTIF directory
    for thestring in (flankbedfile, meme, summit_meme, meme_summary, summit_meme_summary):
        dx_container.move(
            destination=prefolder,
            objects=[thestring['$dnanexus_link']]
        )

#BAM2GFF
    m_downstream = output_map['m_downstream']
    m_upstream = output_map['m_upstream']
    m_genebody = output_map['m_genebody']
    m_promoters = output_map['m_promoters']
    pdf_gene = output_map['pdf_gene']
    pdf_h_gene = output_map['pdf_h_gene']
    png_h_gene = output_map['png_h_gene']
    pdf_promoters = output_map['pdf_promoters']
    pdf_h_promoters = output_map['pdf_h_promoters']
    png_h_promoters = output_map['png_h_promoters']

    prefolder = folder_location + "/" + sample_name + "/BAMDensity_files"
    folder = prefolder + "/matrix"
    dx_container.new_folder(folder, parents=True)
    for thestring in (m_downstream, m_upstream, m_genebody, m_promoters):
        if thestring is not None:
            dx_container.move(
                destination=folder,
                objects=[thestring['$dnanexus_link']]
            )
    for thestring in (pdf_gene, pdf_h_gene, png_h_gene, pdf_promoters,
                      pdf_h_promoters, png_h_promoters):
        if thestring is not None:
            dx_container.move(
                destination=prefolder,
                objects=[thestring['$dnanexus_link']]
            )

#VISUALIZATION
    bigwig = output_map['bigwig']
    norm_wig = output_map['norm_wig']
    tdffile = output_map['tdffile']
    n_bigwig = output_map['n_bigwig']
    n_norm_wig = output_map['n_norm_wig']
    n_tdffile = output_map['n_tdffile']
    a_bigwig = output_map['a_bigwig']
    a_norm_wig = output_map['a_norm_wig']
    a_tdffile = output_map['a_tdffile']

    folder = folder_location + "/" + sample_name + "/PEAKDisplay_files"
    dx_container.new_folder(folder, parents=True)
    for thestring in (bigwig, norm_wig, tdffile,
                      a_bigwig, a_norm_wig, a_tdffile,
                      n_bigwig, n_norm_wig, n_tdffile):
        if thestring is not None:
            dx_container.move(
                destination=folder,
                objects=[thestring['$dnanexus_link']]
            )

    # remove the temporary folder created
    dx_container.remove_folder(folder=temp_folder, recurse=True)
