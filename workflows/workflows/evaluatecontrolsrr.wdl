version 1.0
import "../../seaseq-control.wdl" as sc
import "../../peaseq-control.wdl" as pc
import "../../workflows/tasks/sratoolkit.wdl" as sra

workflow evaluatesrr {
    input {
        File reference
        File? blacklist
        File gtf
        Array[File]? bowtie_index
        Array[File]? motif_databases
        Array[String]? sample_sraid
        Array[File]? sample_R1_fastq
        Array[File]? sample_R2_fastq
        Array[String]? control_sraid
        Array[File]? control_R1_fastq
        Array[File]? control_R2_fastq
        String? results_name
        Boolean run_motifs=true
        Int? insertsize = 600
        String? strandedness = "fr"
        String? results_name
        Boolean run_motifs=true
    }
    
    Boolean paired=true
    if ( defined(sample_sraid) ) {
        Array[String] string_sra = [1] #buffer to allow for sra_id optionality
        Array[String] s_sraid = select_first([sample_sraid, string_sra])
        scatter (eachsra in s_sraid) {
            call sra.srameta {
                input :
                    sra_id=eachsra
            }
            if (! srameta.paired_end){
                Boolean? paired_sample=srameta.paired_end
            }
        } # end scatter each sra
        Boolean paired_sample_m = select_first([paired_sample[0],paired])
    }

    if ( defined(control_sraid) ) {
        Array[String] c_sra = [1] #buffer to allow for sra_id optionality
        Array[String] c_sraid = select_first([control_sraid, c_sra])
        scatter (eachsra in c_sraid) {
            call sra.srameta as c_srameta {
                input :
                    sra_id=eachsra
            }
            if (! c_srameta.paired_end ){
                Boolean? paired_control=c_srameta.paired_end
            }
        } # end scatter each sra
        Boolean paired_control_m = select_first([paired_control[0],paired])
    } # end if control_sra_id

    if (! (paired_control_m && paired_sample_m) ) {
        call sc.seaseq as sc {
            input :
                reference=reference,
                blacklist=blacklist,
                gtf=gtf,
                bowtie_index=bowtie_index,
                motif_databases=motif_databases,
                sample_fastq=sample_R1_fastq,
                control_fastq=control_R1_fastq,
                sample_sraid=sample_sraid,
                control_sraid=control_sraid,
                results_name=results_name,
                run_motifs=run_motifs
        }
    }
    if (paired_sample_m && paired_control_m) {
        call pc.peaseq as pc {
            input :
                reference=reference,
                blacklist=blacklist,
                gtf=gtf,
                bowtie_index=bowtie_index,
                motif_databases=motif_databases,
                sample_R1_fastq=sample_R1_fastq,
                sample_R2_fastq=sample_R2_fastq,
                control_R1_fastq=control_R1_fastq,
                control_R2_fastq=control_R2_fastq,
                sample_sraid=sample_sraid,
                control_sraid=control_sraid,
                insertsize=insertsize,
                strandedness=strandedness,
                results_name=results_name,
                run_motifs=run_motifs
        }
    }

    # Processing OUTPUTs

    output {
        Array[File?]? indv_s_htmlfile = select_first([sc.indv_s_htmlfile, pc.indv_s_htmlfile])
        Array[File?]? indv_s_zipfile = select_first([sc.indv_s_zipfile, pc.indv_s_zipfile])
        Array[File?]? indv_s_bam_htmlfile = select_first([sc.indv_s_bam_htmlfile, pc.indv_s_bam_htmlfile])
        Array[File?]? indv_s_bam_zipfile = select_first([sc.indv_s_bam_zipfile, pc.indv_s_bam_zipfile])
        Array[File?]? indv_c_htmlfile = select_first([sc.indv_c_htmlfile, pc.indv_c_htmlfile])
        Array[File?]? indv_c_zipfile = select_first([sc.indv_c_zipfile, pc.indv_c_zipfile])
        Array[File?]? indv_c_bam_htmlfile = select_first([sc.indv_c_bam_htmlfile, pc.indv_c_bam_htmlfile])
        Array[File?]? indv_c_bam_zipfile = select_first([sc.indv_c_bam_zipfile, pc.indv_c_bam_zipfile])
        File? s_mergebam_htmlfile = select_first([sc.s_mergebam_htmlfile, pc.s_mergebam_htmlfile])
        File? s_mergebam_zipfile = select_first([sc.s_mergebam_zipfile, pc.s_mergebam_zipfile])
        File? c_mergebam_htmlfile = select_first([sc.c_mergebam_htmlfile, pc.c_mergebam_htmlfile])
        File? c_mergebam_zipfile = select_first([sc.c_mergebam_zipfile, pc.c_mergebam_zipfile])
        Array[File?]? indv_sp_bam_htmlfile = pc.indv_sp_bam_htmlfile
        Array[File?]? indv_sp_bam_zipfile = pc.indv_sp_bam_zipfile
        File? sp_mergebam_htmlfile = pc.sp_mergebam_htmlfile
        File? sp_mergebam_zipfile = pc.sp_mergebam_zipfile
        Array[File?]? indv_cp_bam_htmlfile = pc.indv_cp_bam_htmlfile
        Array[File?]? indv_cp_bam_zipfile = pc.indv_cp_bam_zipfile
        File? cp_mergebam_htmlfile = pc.cp_mergebam_htmlfile
        File? cp_mergebam_zipfile = pc.cp_mergebam_zipfile
        File? uno_s_htmlfile = sc.uno_s_htmlfile 
        File? uno_s_zipfile = sc.uno_s_zipfile
        File? uno_s_bam_htmlfile = select_first([sc.uno_s_bam_htmlfile, pc.uno_s_bam_htmlfile])
        File? uno_s_bam_zipfile = select_first([sc.uno_s_bam_zipfile, pc.uno_s_bam_zipfile])
        File? uno_c_htmlfile = sc.uno_c_htmlfile
        File? uno_c_zipfile = sc.uno_c_zipfile
        File? uno_c_bam_htmlfile = select_first([sc.uno_c_bam_htmlfile, pc.uno_c_bam_htmlfile])
        File? uno_c_bam_zipfile = select_first([sc.uno_c_bam_zipfile, pc.uno_c_bam_zipfile])
        Array[File?]? s_metrics_out = select_first([sc.s_metrics_out, pc.s_metrics_out])
        File? uno_s_metrics_out = sc.uno_s_metrics_out
        Array[File?]? c_metrics_out = select_first([sc.c_metrics_out, pc.c_metrics_out])
        File? uno_c_metrics_out = sc.uno_c_metrics_out
        Array[File?]? indv_s_sortedbam = select_first([sc.indv_s_sortedbam, pc.indv_s_sortedbam])
        Array[File?]? indv_s_indexbam = select_first([sc.indv_s_indexbam, pc.indv_s_indexbam])
        Array[File?]? indv_s_bkbam = select_first([sc.indv_s_bkbam, pc.indv_s_bkbam])
        Array[File?]? indv_s_bkindexbam = select_first([sc.indv_s_bkindexbam, pc.indv_s_bkindexbam])
        Array[File?]? indv_s_rmbam = select_first([sc.indv_s_rmbam, pc.indv_s_rmbam])
        Array[File?]? indv_s_rmindexbam = select_first([sc.indv_s_rmindexbam, pc.indv_s_rmindexbam])
        Array[File?]? indv_c_sortedbam = select_first([sc.indv_c_sortedbam, pc.indv_c_sortedbam])
        Array[File?]? indv_c_indexbam = select_first([sc.indv_c_indexbam, pc.indv_c_indexbam])
        Array[File?]? indv_c_bkbam = select_first([sc.indv_c_bkbam, pc.indv_c_bkbam])
        Array[File?]? indv_c_bkindexbam = select_first([sc.indv_c_bkindexbam, pc.indv_c_bkindexbam])
        Array[File?]? indv_c_rmbam = select_first([sc.indv_c_rmbam, pc.indv_c_rmbam])
        Array[File?]? indv_c_rmindexbam = select_first([sc.indv_c_rmindexbam, pc.indv_c_rmindexbam])
        Array[File?]? indv_sp_sortedbam = pc.indv_sp_sortedbam
        Array[File?]? indv_sp_indexbam = pc.indv_sp_indexbam
        Array[File?]? indv_sp_bkbam = pc.indv_sp_bkbam
        Array[File?]? indv_sp_bkindexbam = pc.indv_sp_bkindexbam
        Array[File?]? indv_sp_rmbam = pc.indv_sp_rmbam
        Array[File?]? indv_sp_rmindexbam = pc.indv_sp_rmindexbam
        Array[File?]? indv_cp_sortedbam = pc.indv_cp_sortedbam
        Array[File?]? indv_cp_indexbam = pc.indv_cp_indexbam
        Array[File?]? indv_cp_bkbam = pc.indv_cp_bkbam
        Array[File?]? indv_cp_bkindexbam = pc.indv_cp_bkindexbam
        Array[File?]? indv_cp_rmbam = pc.indv_cp_rmbam
        Array[File?]? indv_cp_rmindexbam = pc.indv_cp_rmindexbam
        File? uno_s_sortedbam = select_first([sc.uno_s_sortedbam, pc.uno_s_sortedbam])
        File? uno_s_indexstatsbam = select_first([sc.uno_s_indexstatsbam, pc.uno_s_indexstatsbam])
        File? uno_s_bkbam = select_first([sc.uno_s_bkbam, pc.uno_s_bkbam])
        File? uno_s_bkindexbam = select_first([sc.uno_s_bkindexbam, pc.uno_s_bkindexbam])
        File? uno_s_rmbam = select_first([sc.uno_s_rmbam, pc.uno_s_rmbam])
        File? uno_s_rmindexbam = select_first([sc.uno_s_rmindexbam, pc.uno_s_rmindexbam])
        File? uno_c_sortedbam = select_first([sc.uno_c_sortedbam, pc.uno_c_sortedbam])
        File? uno_c_indexstatsbam = select_first([sc.uno_c_indexstatsbam, pc.uno_c_indexstatsbam])
        File? uno_c_bkbam = select_first([sc.uno_c_bkbam, pc.uno_c_bkbam])
        File? uno_c_bkindexbam = select_first([sc.uno_c_bkindexbam, pc.uno_c_bkindexbam])
        File? uno_c_rmbam = select_first([sc.uno_c_rmbam, pc.uno_c_rmbam])
        File? uno_c_rmindexbam = select_first([sc.uno_c_rmindexbam, pc.uno_c_rmindexbam])
        File? s_mergebamfile = select_first([sc.s_mergebamfile, pc.s_mergebamfile])
        File? s_mergebamindex = select_first([sc.s_mergebamindex, pc.s_mergebamindex])
        File? s_bkbam = select_first([sc.s_bkbam, pc.s_bkbam])
        File? s_bkindexbam = select_first([sc.s_bkindexbam, pc.s_bkindexbam])
        File? s_rmbam = select_first([sc.s_rmbam, pc.s_rmbam])
        File? s_rmindexbam = select_first([sc.s_rmindexbam, pc.s_rmindexbam])
        File? c_mergebamfile = select_first([sc.c_mergebamfile, pc.c_mergebamfile])
        File? c_mergebamindex = select_first([sc.c_mergebamindex, pc.c_mergebamindex])
        File? c_bkbam = select_first([sc.c_bkbam, pc.c_bkbam])
        File? c_bkindexbam = select_first([sc.c_bkindexbam, pc.c_bkindexbam])
        File? c_rmbam = select_first([sc.c_rmbam, pc.c_rmbam])
        File? c_rmindexbam = select_first([sc.c_rmindexbam, pc.c_rmindexbam])
        File? sp_mergebamfile = pc.sp_mergebamfile
        File? sp_mergebamindex = pc.sp_mergebamindex
        File? sp_bkbam = pc.sp_bkbam
        File? sp_bkindexbam = pc.sp_bkindexbam
        File? sp_rmbam = pc.sp_rmbam
        File? sp_rmindexbam = pc.sp_rmindexbam
        File? cp_mergebamfile = pc.cp_mergebamfile
        File? cp_mergebamindex = pc.cp_mergebamindex
        File? cp_bkbam = pc.cp_bkbam
        File? cp_bkindexbam = pc.cp_bkindexbam
        File? cp_rmbam = pc.cp_rmbam
        File? cp_rmindexbam = pc.cp_rmindexbam
        File? s_fragments_bam = pc.s_fragments_bam
        File? s_fragments_indexbam = pc.s_fragments_indexbam
        File? c_fragments_bam = pc.c_fragments_bam
        File? c_fragments_indexbam = pc.c_fragments_indexbam
        File? peakbedfile = select_first([sc.peakbedfile, pc.peakbedfile])
        File? peakxlsfile = select_first([sc.peakxlsfile, pc.peakxlsfile])
        File? negativexlsfile = select_first([sc.negativexlsfile, pc.negativexlsfile])
        File? summitsfile = select_first([sc.summitsfile, pc.summitsfile])
        File? wigfile = select_first([sc.wigfile, pc.wigfile])
        File? ctrlwigfile = select_first([sc.ctrlwigfile, pc.ctrlwigfile])
        File? all_peakbedfile = select_first([sc.all_peakbedfile, pc.all_peakbedfile])
        File? all_peakxlsfile = select_first([sc.all_peakxlsfile, pc.all_peakxlsfile])
        File? all_negativexlsfile = select_first([sc.all_negativexlsfile, pc.all_negativexlsfile])
        File? all_summitsfile = select_first([sc.all_summitsfile, pc.all_summitsfile])
        File? all_wigfile = select_first([sc.all_wigfile, pc.all_wigfile])
        File? all_ctrlwigfile = select_first([sc.all_ctrlwigfile, pc.all_ctrlwigfile])
        File? nm_peakbedfile = select_first([sc.nm_peakbedfile, pc.nm_peakbedfile])
        File? nm_peakxlsfile = select_first([sc.nm_peakxlsfile, pc.nm_peakxlsfile])
        File? nm_negativexlsfile = select_first([sc.nm_negativexlsfile, pc.nm_negativexlsfile])
        File? nm_summitsfile = select_first([sc.nm_summitsfile, pc.nm_summitsfile])
        File? nm_wigfile = select_first([sc.nm_wigfile, pc.nm_wigfile])
        File? nm_ctrlwigfile = select_first([sc.nm_ctrlwigfile, pc.nm_ctrlwigfile])
        File? readme_peaks = select_first([sc.readme_peaks, pc.readme_peaks])
        File? only_s_peakbedfile = select_first([sc.only_s_peakbedfile, pc.only_s_peakbedfile])
        File? only_s_peakxlsfile = select_first([sc.only_s_peakxlsfile, pc.only_s_peakxlsfile])
        File? only_s_summitsfile = select_first([sc.only_s_summitsfile, pc.only_s_summitsfile])
        File? only_s_wigfile = select_first([sc.only_s_wigfile, pc.only_s_wigfile])
        File? only_c_peakbedfile = select_first([sc.only_c_peakbedfile, pc.only_c_peakbedfile])
        File? only_c_peakxlsfile = select_first([sc.only_c_peakxlsfile, pc.only_c_peakxlsfile])
        File? only_c_summitsfile = select_first([sc.only_c_summitsfile, pc.only_c_summitsfile])
        File? only_c_wigfile = select_first([sc.only_c_wigfile, pc.only_c_wigfile])
        File? only_sp_peakbedfile = pc.only_sp_peakbedfile
        File? only_sp_peakxlsfile = pc.only_sp_peakxlsfile
        File? only_sp_summitsfile = pc.only_sp_summitsfile
        File? only_sp_wigfile = pc.only_sp_wigfile
        File? only_cp_peakbedfile = pc.only_cp_peakbedfile
        File? only_cp_peakxlsfile = pc.only_cp_peakxlsfile
        File? only_cp_summitsfile = pc.only_cp_summitsfile
        File? only_cp_wigfile = pc.only_cp_wigfile
        File? sp_peakbedfile = pc.sp_peakbedfile
        File? sp_peakxlsfile = pc.sp_peakxlsfile
        File? sp_negativexlsfile = pc.sp_negativexlsfile
        File? sp_summitsfile = pc.sp_summitsfile
        File? sp_wigfile = pc.sp_wigfile
        File? sp_ctrlwigfile = pc.sp_ctrlwigfile
        File? sp_all_peakbedfile = pc.sp_all_peakbedfile
        File? sp_all_peakxlsfile = pc.sp_all_peakxlsfile
        File? sp_all_negativexlsfile = pc.sp_all_negativexlsfile
        File? sp_all_summitsfile = pc.sp_all_summitsfile
        File? sp_all_wigfile = pc.sp_all_wigfile
        File? sp_all_ctrlwigfile = pc.sp_all_ctrlwigfile
        File? sp_nm_peakbedfile = pc.sp_nm_peakbedfile
        File? sp_nm_peakxlsfile = pc.sp_nm_peakxlsfile
        File? sp_nm_negativexlsfile = pc.sp_nm_negativexlsfile
        File? sp_nm_summitsfile = pc.sp_nm_summitsfile
        File? sp_nm_wigfile = pc.sp_nm_wigfile
        File? sp_nm_ctrlwigfile = pc.sp_nm_ctrlwigfile
        File? sp_readme_peaks = pc.sp_readme_peaks
        File? scoreisland = select_first([sc.scoreisland, pc.scoreisland])
        File? sicer_wigfile = select_first([sc.sicer_wigfile, pc.sicer_wigfile])
        File? sicer_summary = select_first([sc.sicer_summary, pc.sicer_summary])
        File? sicer_fdrisland = select_first([sc.sicer_fdrisland, pc.sicer_fdrisland])
        File? sp_scoreisland = pc.sp_scoreisland
        File? sp_sicer_wigfile = pc.sp_sicer_wigfile
        File? sp_sicer_summary = pc.sp_sicer_summary
        File? sp_sicer_fdrisland = pc.sp_sicer_fdrisland
        File? pngfile = select_first([sc.pngfile, pc.pngfile])
        File? mapped_union = select_first([sc.mapped_union, pc.mapped_union])
        File? mapped_stitch = select_first([sc.mapped_stitch, pc.mapped_stitch])
        File? enhancers = select_first([sc.enhancers, pc.enhancers])
        File? super_enhancers = select_first([sc.super_enhancers, pc.super_enhancers])
        File? gff_file = select_first([sc.gff_file, pc.gff_file])
        File? gff_union = select_first([sc.gff_union, pc.gff_union])
        File? union_enhancers = select_first([sc.union_enhancers, pc.union_enhancers])
        File? stitch_enhancers = select_first([sc.stitch_enhancers, pc.stitch_enhancers])
        File? e_to_g_enhancers = select_first([sc.e_to_g_enhancers, pc.e_to_g_enhancers])
        File? g_to_e_enhancers = select_first([sc.g_to_e_enhancers, pc.g_to_e_enhancers])
        File? e_to_g_super_enhancers = select_first([sc.e_to_g_super_enhancers, pc.e_to_g_super_enhancers])
        File? g_to_e_super_enhancers = select_first([sc.g_to_e_super_enhancers, pc.g_to_e_super_enhancers])
        File? sp_pngfile = pc.sp_pngfile
        File? sp_mapped_union = pc.sp_mapped_union
        File? sp_mapped_stitch = pc.sp_mapped_stitch
        File? sp_enhancers = pc.sp_enhancers
        File? sp_super_enhancers = pc.sp_super_enhancers
        File? sp_gff_file = pc.sp_gff_file
        File? sp_gff_union = pc.sp_gff_union
        File? sp_union_enhancers = pc.sp_union_enhancers
        File? sp_stitch_enhancers = pc.sp_stitch_enhancers
        File? sp_e_to_g_enhancers = pc.sp_e_to_g_enhancers
        File? sp_g_to_e_enhancers = pc.sp_g_to_e_enhancers
        File? sp_e_to_g_super_enhancers = pc.sp_e_to_g_super_enhancers
        File? sp_g_to_e_super_enhancers = pc.sp_g_to_e_super_enhancers
        File? flankbedfile = select_first([sc.flankbedfile, pc.flankbedfile])
        File? ame_tsv = select_first([sc.ame_tsv, pc.ame_tsv])
        File? ame_html = select_first([sc.ame_html, pc.ame_html])
        File? ame_seq = select_first([sc.ame_seq, pc.ame_seq])
        File? meme = select_first([sc.meme, pc.meme])
        File? meme_summary = select_first([sc.meme_summary, pc.meme_summary])
        File? summit_ame_tsv = select_first([sc.summit_ame_tsv, pc.summit_ame_tsv])
        File? summit_ame_html = select_first([sc.summit_ame_html, pc.summit_ame_html])
        File? summit_ame_seq = select_first([sc.summit_ame_seq, pc.summit_ame_seq])
        File? summit_meme = select_first([sc.summit_meme, pc.summit_meme])
        File? summit_meme_summary = select_first([sc.summit_meme_summary, pc.summit_meme_summary])
        File? sp_flankbedfile = pc.sp_flankbedfile
        File? sp_ame_tsv = pc.sp_ame_tsv
        File? sp_ame_html = pc.sp_ame_html
        File? sp_ame_seq = pc.sp_ame_seq
        File? sp_meme = pc.sp_meme
        File? sp_meme_summary = pc.sp_meme_summary
        File? sp_summit_ame_tsv = pc.sp_summit_ame_tsv
        File? sp_summit_ame_html = pc.sp_summit_ame_html
        File? sp_summit_ame_seq = pc.sp_summit_ame_seq
        File? sp_summit_meme = pc.sp_summit_meme
        File? sp_summit_meme_summary = pc.sp_summit_meme_summary
        File? s_matrices = select_first([sc.s_matrices, pc.s_matrices])
        File? c_matrices = select_first([sc.c_matrices, pc.c_matrices])
        File? densityplot = select_first([sc.densityplot, pc.densityplot])
        File? pdf_gene = select_first([sc.pdf_gene, pc.pdf_gene])
        File? pdf_h_gene = select_first([sc.pdf_h_gene, pc.pdf_h_gene])
        File? png_h_gene = select_first([sc.png_h_gene, pc.png_h_gene])
        File? jpg_h_gene = select_first([sc.jpg_h_gene, pc.jpg_h_gene])
        File? pdf_promoters = select_first([sc.pdf_promoters, pc.pdf_promoters])
        File? pdf_h_promoters = select_first([sc.pdf_h_promoters, pc.pdf_h_promoters])
        File? png_h_promoters = select_first([sc.png_h_promoters, pc.png_h_promoters])
        File? jpg_h_promoters = select_first([sc.jpg_h_promoters, pc.jpg_h_promoters])
        File? sp_s_matrices = pc.sp_s_matrices
        File? sp_c_matrices = pc.sp_c_matrices
        File? sp_densityplot = pc.sp_densityplot
        File? sp_pdf_gene = pc.sp_pdf_gene
        File? sp_pdf_h_gene = pc.sp_pdf_h_gene
        File? sp_png_h_gene = pc.sp_png_h_gene
        File? sp_jpg_h_gene = pc.sp_jpg_h_gene
        File? sp_pdf_promoters = pc.sp_pdf_promoters
        File? sp_pdf_h_promoters = pc.sp_pdf_h_promoters
        File? sp_png_h_promoters = pc.sp_png_h_promoters
        File? sp_jpg_h_promoters = pc.sp_jpg_h_promoters
        File? peak_promoters = select_first([sc.peak_promoters, pc.peak_promoters])
        File? peak_genebody = select_first([sc.peak_genebody, pc.peak_genebody])
        File? peak_window = select_first([sc.peak_window, pc.peak_window])
        File? peak_closest = select_first([sc.peak_closest, pc.peak_closest])
        File? peak_comparison = select_first([sc.peak_comparison, pc.peak_comparison])
        File? gene_comparison = select_first([sc.gene_comparison, pc.gene_comparison])
        File? pdf_comparison = select_first([sc.pdf_comparison, pc.pdf_comparison])
        File? all_peak_promoters = select_first([sc.all_peak_promoters, pc.all_peak_promoters])
        File? all_peak_genebody = select_first([sc.all_peak_genebody, pc.all_peak_genebody])
        File? all_peak_window = select_first([sc.all_peak_window, pc.all_peak_window])
        File? all_peak_closest = select_first([sc.all_peak_closest, pc.all_peak_closest])
        File? all_peak_comparison = select_first([sc.all_peak_comparison, pc.all_peak_comparison])
        File? all_gene_comparison = select_first([sc.all_gene_comparison, pc.all_gene_comparison])
        File? all_pdf_comparison = select_first([sc.all_pdf_comparison, pc.all_pdf_comparison])
        File? nomodel_peak_promoters = select_first([sc.nomodel_peak_promoters, pc.nomodel_peak_promoters])
        File? nomodel_peak_genebody = select_first([sc.nomodel_peak_genebody, pc.nomodel_peak_genebody])
        File? nomodel_peak_window = select_first([sc.nomodel_peak_window, pc.nomodel_peak_window])
        File? nomodel_peak_closest = select_first([sc.nomodel_peak_closest, pc.nomodel_peak_closest])
        File? nomodel_peak_comparison = select_first([sc.nomodel_peak_comparison, pc.nomodel_peak_comparison])
        File? nomodel_gene_comparison = select_first([sc.nomodel_gene_comparison, pc.nomodel_gene_comparison])
        File? nomodel_pdf_comparison = select_first([sc.nomodel_pdf_comparison, pc.nomodel_pdf_comparison])
        File? sicer_peak_promoters = select_first([sc.sicer_peak_promoters, pc.sicer_peak_promoters])
        File? sicer_peak_genebody = select_first([sc.sicer_peak_genebody, pc.sicer_peak_genebody])
        File? sicer_peak_window = select_first([sc.sicer_peak_window, pc.sicer_peak_window])
        File? sicer_peak_closest = select_first([sc.sicer_peak_closest, pc.sicer_peak_closest])
        File? sicer_peak_comparison = select_first([sc.sicer_peak_comparison, pc.sicer_peak_comparison])
        File? sicer_gene_comparison = select_first([sc.sicer_gene_comparison, pc.sicer_gene_comparison])
        File? sicer_pdf_comparison = select_first([sc.sicer_pdf_comparison, pc.sicer_pdf_comparison])
        File? sp_peak_promoters = pc.sp_peak_promoters
        File? sp_peak_genebody =  pc.sp_peak_genebody
        File? sp_peak_window = pc.sp_peak_window
        File? sp_peak_closest = pc.sp_peak_closest
        File? sp_peak_comparison = pc.sp_peak_comparison
        File? sp_gene_comparison = pc.sp_gene_comparison
        File? sp_pdf_comparison = pc.sp_pdf_comparison
        File? sp_all_peak_promoters = pc.sp_all_peak_promoters
        File? sp_all_peak_genebody =  pc.sp_all_peak_genebody
        File? sp_all_peak_window = pc.sp_all_peak_window
        File? sp_all_peak_closest = pc.sp_all_peak_closest
        File? sp_all_peak_comparison = pc.sp_all_peak_comparison
        File? sp_all_gene_comparison = pc.sp_all_gene_comparison
        File? sp_all_pdf_comparison = pc.sp_all_pdf_comparison
        File? sp_nomodel_peak_promoters = pc.sp_nomodel_peak_promoters
        File? sp_nomodel_peak_genebody =  pc.sp_nomodel_peak_genebody
        File? sp_nomodel_peak_window = pc.sp_nomodel_peak_window
        File? sp_nomodel_peak_closest = pc.sp_nomodel_peak_closest
        File? sp_nomodel_peak_comparison = pc.sp_nomodel_peak_comparison
        File? sp_nomodel_gene_comparison = pc.sp_nomodel_gene_comparison
        File? sp_nomodel_pdf_comparison = pc.sp_nomodel_pdf_comparison
        File? sp_sicer_peak_promoters = pc.sp_sicer_peak_promoters
        File? sp_sicer_peak_genebody =  pc.sp_sicer_peak_genebody
        File? sp_sicer_peak_window = pc.sp_sicer_peak_window
        File? sp_sicer_peak_closest = pc.sp_sicer_peak_closest
        File? sp_sicer_peak_comparison = pc.sp_sicer_peak_comparison
        File? sp_sicer_gene_comparison = pc.sp_sicer_gene_comparison
        File? sp_sicer_pdf_comparison = pc.sp_sicer_pdf_comparison
        File? bigwig = select_first([sc.bigwig, pc.bigwig])
        File? norm_wig = select_first([sc.norm_wig, pc.norm_wig])
        File? tdffile = select_first([sc.tdffile, pc.tdffile])
        File? n_bigwig = select_first([sc.n_bigwig, pc.n_bigwig])
        File? n_norm_wig = select_first([sc.n_norm_wig, pc.n_norm_wig])
        File? n_tdffile = select_first([sc.n_tdffile, pc.n_tdffile])
        File? a_bigwig = select_first([sc.a_bigwig, pc.a_bigwig])
        File? a_norm_wig = select_first([sc.a_norm_wig, pc.a_norm_wig])
        File? a_tdffile = select_first([sc.a_tdffile, pc.a_tdffile])
        File? c_bigwig = select_first([sc.c_bigwig, pc.c_bigwig])
        File? c_norm_wig = select_first([sc.c_norm_wig, pc.c_norm_wig])
        File? c_tdffile = select_first([sc.c_tdffile, pc.c_tdffile])
        File? c_n_bigwig = select_first([sc.c_n_bigwig, pc.c_n_bigwig])
        File? c_n_norm_wig = select_first([sc.c_n_norm_wig, pc.c_n_norm_wig])
        File? c_n_tdffile = select_first([sc.c_n_tdffile, pc.c_n_tdffile])
        File? c_a_bigwig = select_first([sc.c_a_bigwig, pc.c_a_bigwig])
        File? c_a_norm_wig = select_first([sc.c_a_norm_wig, pc.c_a_norm_wig])
        File? c_a_tdffile = select_first([sc.c_a_tdffile, pc.c_a_tdffile])
        File? s_bigwig = select_first([sc.s_bigwig, pc.s_bigwig])
        File? s_norm_wig = select_first([sc.s_norm_wig, pc.s_norm_wig])
        File? s_tdffile = select_first([sc.s_tdffile, pc.s_tdffile])
        File? sp_bigwig = pc.sp_bigwig
        File? sp_norm_wig = pc.sp_norm_wig
        File? sp_tdffile = pc.sp_tdffile
        File? sp_n_bigwig = pc.sp_n_bigwig
        File? sp_n_norm_wig = pc.sp_n_norm_wig
        File? sp_n_tdffile = pc.sp_n_tdffile
        File? sp_a_bigwig = pc.sp_a_bigwig
        File? sp_a_norm_wig = pc.sp_a_norm_wig
        File? sp_a_tdffile = pc.sp_a_tdffile
        File? cp_bigwig = pc.cp_bigwig
        File? cp_norm_wig = pc.cp_norm_wig
        File? cp_tdffile = pc.cp_tdffile
        File? cp_n_bigwig = pc.cp_n_bigwig
        File? cp_n_norm_wig = pc.cp_n_norm_wig
        File? cp_n_tdffile = pc.cp_n_tdffile
        File? cp_a_bigwig = pc.cp_a_bigwig
        File? cp_a_norm_wig = pc.cp_a_norm_wig
        File? cp_a_tdffile = pc.cp_a_tdffile
        File? sp_s_bigwig = pc.sp_s_bigwig
        File? sp_s_norm_wig = pc.sp_s_norm_wig
        File? sp_s_tdffile = pc.sp_s_tdffile
        File? sf_bigwig = pc.sf_bigwig
        File? sf_tdffile = pc.sf_tdffile
        File? sf_wigfile = pc.sf_wigfile
        File? cf_bigwig = pc.cf_bigwig
        File? cf_tdffile = pc.cf_tdffile
        File? cf_wigfile = pc.cf_wigfile
        Array[File?]? s_qc_statsfile = select_first([sc.s_qc_statsfile, pc.s_qc_statsfile])
        Array[File?]? s_qc_htmlfile = select_first([sc.s_qc_htmlfile, pc.s_qc_htmlfile])
        Array[File?]? s_qc_textfile = select_first([sc.s_qc_textfile, pc.s_qc_textfile])
        File? s_qc_mergehtml = select_first([sc.s_qc_mergehtml, pc.s_qc_mergehtml])
        Array[File?]? c_qc_statsfile = select_first([sc.c_qc_statsfile, pc.c_qc_statsfile])
        Array[File?]? c_qc_htmlfile = select_first([sc.c_qc_htmlfile, pc.c_qc_htmlfile])
        Array[File?]? c_qc_textfile = select_first([sc.c_qc_textfile, pc.c_qc_textfile])
        File? c_qc_mergehtml = select_first([sc.c_qc_mergehtml, pc.c_qc_mergehtml])
        Array[File?]? sp_qc_statsfile = pc.sp_qc_statsfile
        Array[File?]? sp_qc_htmlfile = pc.sp_qc_htmlfile
        Array[File?]? sp_qc_textfile = pc.sp_qc_textfile
        Array[File?]? cp_qc_statsfile = pc.cp_qc_statsfile
        Array[File?]? cp_qc_htmlfile = pc.cp_qc_htmlfile
        Array[File?]? cp_qc_textfile = pc.cp_qc_textfile
        File? s_uno_statsfile = select_first([sc.s_uno_statsfile, pc.s_uno_statsfile])
        File? s_uno_htmlfile = select_first([sc.s_uno_htmlfile, pc.s_uno_htmlfile])
        File? s_uno_textfile = select_first([sc.s_uno_textfile, pc.s_uno_textfile])
        File? c_uno_statsfile = select_first([sc.c_uno_statsfile, pc.c_uno_statsfile])
        File? c_uno_htmlfile = select_first([sc.c_uno_htmlfile, pc.c_uno_htmlfile])
        File? c_uno_textfile = select_first([sc.c_uno_textfile, pc.c_uno_textfile])
        File? statsfile = select_first([sc.statsfile, pc.statsfile])
        File? htmlfile = select_first([sc.htmlfile, pc.htmlfile])
        File? textfile = select_first([sc.textfile, pc.textfile])
        File? s_statsfile = pc.s_statsfile
        File? s_htmlfile = pc.s_htmlfile
        File? s_textfile = pc.s_textfile
        File? c_statsfile = pc.c_statsfile
        File? c_htmlfile = pc.c_htmlfile
        File? c_textfile = pc.c_textfile
        File? sp_statsfile = pc.sp_statsfile
        File? sp_htmlfile = pc.sp_htmlfile
        File? sp_textfile = pc.sp_textfile
        File? s_summaryhtml = pc.s_summaryhtml
        File? s_summarystats = pc.s_summarystats
        File? s_summarytxt = pc.s_summarytxt
        File? c_summaryhtml = pc.c_summaryhtml
        File? c_summarystats = pc.c_summarystats
        File? c_summarytxt = pc.c_summarytxt
        File? s_fragsize = pc.s_fragsize
        File? c_fragsize = pc.c_fragsize
        File? summaryhtml = select_first([sc.summaryhtml, pc.summaryhtml])
        File? summarytxt = select_first([sc.summarytxt, pc.summarytxt])
        File? sp_summaryhtml = pc.sp_summaryhtml
        File? sp_summarytxt = pc.sp_summarytxt
    }
}
