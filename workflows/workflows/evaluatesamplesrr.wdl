version 1.0
import "../../seaseq-case.wdl" as ss
import "../../peaseq-case.wdl" as ps
import "../../workflows/tasks/sratoolkit.wdl" as sra

workflow evaluatesrr {
    input {
        File reference
        File? spikein_reference
        File? blacklist
        File gtf
        Array[File]? bowtie_index
        Array[File]? spikein_bowtie_index
        Array[File]? motif_databases
        Array[String]? sample_sraid
        Array[File]? sample_R1_fastq
        Array[File]? sample_R2_fastq
        String? results_name
        Boolean run_motifs=true
        Int? insertsize = 600
        String? strandedness = "fr"
        Boolean paired=true
    }
    
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

    if (!paired_sample_m) {
        call ss.seaseq as ss {
            input :
                reference=reference,
                spikein_reference=spikein_reference
                blacklist=blacklist,
                gtf=gtf,
                bowtie_index=bowtie_index,
                spikein_bowtie_index=spikein_bowtie_index,
                motif_databases=motif_databases,
                sample_fastq=sample_R1_fastq,
                sample_sraid=sample_sraid,
                results_name=results_name,
                run_motifs=run_motifs
        }
    }

    if (paired_sample_m) {
        call ps.peaseq as ps {
            input :
                reference=reference,
                spikein_reference=spikein_reference
                blacklist=blacklist,
                gtf=gtf,
                bowtie_index=bowtie_index,
                spikein_bowtie_index=spikein_bowtie_index,
                motif_databases=motif_databases,
                sample_R1_fastq=sample_R1_fastq,
                sample_R2_fastq=sample_R2_fastq,
                sample_sraid=sample_sraid,
                insertsize=insertsize,
                strandedness=strandedness,
                results_name=results_name,
                run_motifs=run_motifs
        }
    }

    # Processing OUTPUTs

    output {
        Array[File?]? spikein_indv_s_htmlfile = if paired_sample_m then ps.spikein_indv_s_htmlfile else ss.spikein_indv_s_htmlfile
        Array[File?]? spikein_indv_s_zipfile = if paired_sample_m then ps.spikein_indv_s_zipfile else ss.spikein_indv_s_zipfile
        Array[File?]? spikein_s_metrics_out = if paired_sample_m then ps.spikein_s_metrics_out else ss.spikein_s_metrics_out
                
        Array[File?]? indv_s_htmlfile = if paired_sample_m then ps.indv_s_htmlfile else ss.indv_s_htmlfile
        Array[File?]? indv_s_zipfile = if paired_sample_m then ps.indv_s_zipfile else ss.indv_s_zipfile
        Array[File?]? indv_s_bam_htmlfile = if paired_sample_m then ps.indv_s_bam_htmlfile else ss.indv_s_bam_htmlfile
        Array[File?]? indv_s_bam_zipfile = if paired_sample_m then ps.indv_s_bam_zipfile else ss.indv_s_bam_zipfile
        File? s_mergebam_htmlfile = if paired_sample_m then ps.s_mergebam_htmlfile else ss.s_mergebam_htmlfile
        File? s_mergebam_zipfile = if paired_sample_m then ps.s_mergebam_zipfile else ss.s_mergebam_zipfile
        Array[File?]? indv_sp_bam_htmlfile = ps.indv_sp_bam_htmlfile
        Array[File?]? indv_sp_bam_zipfile = ps.indv_sp_bam_zipfile
        File? sp_mergebam_htmlfile = ps.sp_mergebam_htmlfile
        File? sp_mergebam_zipfile = ps.sp_mergebam_zipfile
        File? uno_s_htmlfile = ss.uno_s_htmlfile
        File? uno_s_zipfile = ss.uno_s_zipfile
        File? uno_s_bam_htmlfile = if paired_sample_m then ps.uno_s_bam_htmlfile else ss.uno_s_bam_htmlfile
        File? uno_s_bam_zipfile = if paired_sample_m then ps.uno_s_bam_zipfile else ss.uno_s_bam_zipfile
        Array[File?]? s_metrics_out = if paired_sample_m then ps.s_metrics_out else ss.s_metrics_out
        File? uno_s_metrics_out = ss.uno_s_metrics_out
        Array[File?]? indv_s_sortedbam = if paired_sample_m then ps.indv_s_sortedbam else ss.indv_s_sortedbam
        Array[File?]? indv_s_indexbam = if paired_sample_m then ps.indv_s_indexbam else ss.indv_s_indexbam
        Array[File?]? indv_s_bkbam = if paired_sample_m then ps.indv_s_bkbam else ss.indv_s_bkbam
        Array[File?]? indv_s_bkindexbam = if paired_sample_m then ps.indv_s_bkindexbam else ss.indv_s_bkindexbam
        Array[File?]? indv_s_rmbam = if paired_sample_m then ps.indv_s_rmbam else ss.indv_s_rmbam
        Array[File?]? indv_s_rmindexbam = if paired_sample_m then ps.indv_s_rmindexbam else ss.indv_s_rmindexbam
        Array[File?]? indv_sp_sortedbam = ps.indv_sp_sortedbam
        Array[File?]? indv_sp_indexbam = ps.indv_sp_indexbam
        Array[File?]? indv_sp_bkbam = ps.indv_sp_bkbam
        Array[File?]? indv_sp_bkindexbam = ps.indv_sp_bkindexbam
        Array[File?]? indv_sp_rmbam = ps.indv_sp_rmbam
        Array[File?]? indv_sp_rmindexbam = ps.indv_sp_rmindexbam
        File? uno_s_sortedbam = if paired_sample_m then ps.uno_s_sortedbam else ss.uno_s_sortedbam
        File? uno_s_indexstatsbam = if paired_sample_m then ps.uno_s_indexstatsbam else ss.uno_s_indexstatsbam
        File? uno_s_bkbam = if paired_sample_m then ps.uno_s_bkbam else ss.uno_s_bkbam
        File? uno_s_bkindexbam = if paired_sample_m then ps.uno_s_bkindexbam else ss.uno_s_bkindexbam
        File? uno_s_rmbam = if paired_sample_m then ps.uno_s_rmbam else ss.uno_s_rmbam
        File? uno_s_rmindexbam = if paired_sample_m then ps.uno_s_rmindexbam else ss.uno_s_rmindexbam
        File? s_mergebamfile = if paired_sample_m then ps.s_mergebamfile else ss.s_mergebamfile
        File? s_mergebamindex = if paired_sample_m then ps.s_mergebamindex else ss.s_mergebamindex
        File? s_bkbam = if paired_sample_m then ps.s_bkbam else ss.s_bkbam
        File? s_bkindexbam = if paired_sample_m then ps.s_bkindexbam else ss.s_bkindexbam
        File? s_rmbam = if paired_sample_m then ps.s_rmbam else ss.s_rmbam
        File? s_rmindexbam = if paired_sample_m then ps.s_rmindexbam else ss.s_rmindexbam
        File? sp_mergebamfile = ps.sp_mergebamfile
        File? sp_mergebamindex = ps.sp_mergebamindex
        File? sp_bkbam = ps.sp_bkbam
        File? sp_bkindexbam = ps.sp_bkindexbam
        File? sp_rmbam = ps.sp_rmbam
        File? sp_rmindexbam = ps.sp_rmindexbam
        File? s_fragments_bam = ps.s_fragments_bam
        File? s_fragments_indexbam = ps.s_fragments_indexbam
        File? peakbedfile = if paired_sample_m then ps.peakbedfile else ss.peakbedfile
        File? peakxlsfile = if paired_sample_m then ps.peakxlsfile else ss.peakxlsfile
        File? negativexlsfile = if paired_sample_m then ps.negativexlsfile else ss.negativexlsfile
        File? summitsfile = if paired_sample_m then ps.summitsfile else ss.summitsfile
        File? wigfile = if paired_sample_m then ps.wigfile else ss.wigfile
        File? all_peakbedfile = if paired_sample_m then ps.all_peakbedfile else ss.all_peakbedfile
        File? all_peakxlsfile = if paired_sample_m then ps.all_peakxlsfile else ss.all_peakxlsfile
        File? all_negativexlsfile = if paired_sample_m then ps.all_negativexlsfile else ss.all_negativexlsfile
        File? all_summitsfile = if paired_sample_m then ps.all_summitsfile else ss.all_summitsfile
        File? all_wigfile = if paired_sample_m then ps.all_wigfile else ss.all_wigfile
        File? nm_peakbedfile = if paired_sample_m then ps.nm_peakbedfile else ss.nm_peakbedfile
        File? nm_peakxlsfile = if paired_sample_m then ps.nm_peakxlsfile else ss.nm_peakxlsfile
        File? nm_negativexlsfile = if paired_sample_m then ps.nm_negativexlsfile else ss.nm_negativexlsfile
        File? nm_summitsfile = if paired_sample_m then ps.nm_summitsfile else ss.nm_summitsfile
        File? nm_wigfile = if paired_sample_m then ps.nm_wigfile else ss.nm_wigfile
        File? readme_peaks = if paired_sample_m then ps.readme_peaks else ss.readme_peaks
        File? sp_peakbedfile = ps.sp_peakbedfile
        File? sp_peakxlsfile = ps.sp_peakxlsfile
        File? sp_negativexlsfile = ps.sp_negativexlsfile
        File? sp_summitsfile = ps.sp_summitsfile
        File? sp_wigfile = ps.sp_wigfile
        File? sp_all_peakbedfile = ps.sp_all_peakbedfile
        File? sp_all_peakxlsfile = ps.sp_all_peakxlsfile
        File? sp_all_negativexlsfile = ps.sp_all_negativexlsfile
        File? sp_all_summitsfile = ps.sp_all_summitsfile
        File? sp_all_wigfile = ps.sp_all_wigfile
        File? sp_nm_peakbedfile = ps.sp_nm_peakbedfile
        File? sp_nm_peakxlsfile = ps.sp_nm_peakxlsfile
        File? sp_nm_negativexlsfile = ps.sp_nm_negativexlsfile
        File? sp_nm_summitsfile = ps.sp_nm_summitsfile
        File? sp_nm_wigfile = ps.sp_nm_wigfile
        File? sp_readme_peaks = ps.sp_readme_peaks
        File? scoreisland = if paired_sample_m then ps.scoreisland else ss.scoreisland
        File? sicer_wigfile = if paired_sample_m then ps.sicer_wigfile else ss.sicer_wigfile
        File? sp_scoreisland = ps.sp_scoreisland
        File? sp_sicer_wigfile = ps.sp_sicer_wigfile
        File? pngfile = if paired_sample_m then ps.pngfile else ss.pngfile
        File? mapped_union = if paired_sample_m then ps.mapped_union else ss.mapped_union
        File? mapped_stitch = if paired_sample_m then ps.mapped_stitch else ss.mapped_stitch
        File? enhancers = if paired_sample_m then ps.enhancers else ss.enhancers
        File? super_enhancers = if paired_sample_m then ps.super_enhancers else ss.super_enhancers
        File? gff_file = if paired_sample_m then ps.gff_file else ss.gff_file
        File? gff_union = if paired_sample_m then ps.gff_union else ss.gff_union
        File? union_enhancers = if paired_sample_m then ps.union_enhancers else ss.union_enhancers
        File? stitch_enhancers = if paired_sample_m then ps.stitch_enhancers else ss.stitch_enhancers
        File? e_to_g_enhancers = if paired_sample_m then ps.e_to_g_enhancers else ss.e_to_g_enhancers
        File? g_to_e_enhancers = if paired_sample_m then ps.g_to_e_enhancers else ss.g_to_e_enhancers
        File? e_to_g_super_enhancers = if paired_sample_m then ps.e_to_g_super_enhancers else ss.e_to_g_super_enhancers
        File? g_to_e_super_enhancers = if paired_sample_m then ps.g_to_e_super_enhancers else ss.g_to_e_super_enhancers
        File? sp_pngfile = ps.sp_pngfile
        File? sp_mapped_union = ps.sp_mapped_union
        File? sp_mapped_stitch = ps.sp_mapped_stitch
        File? sp_enhancers = ps.sp_enhancers
        File? sp_super_enhancers = ps.sp_super_enhancers
        File? sp_gff_file = ps.sp_gff_file
        File? sp_gff_union = ps.sp_gff_union
        File? sp_union_enhancers = ps.sp_union_enhancers
        File? sp_stitch_enhancers = ps.sp_stitch_enhancers
        File? sp_e_to_g_enhancers = ps.sp_e_to_g_enhancers
        File? sp_g_to_e_enhancers = ps.sp_g_to_e_enhancers
        File? sp_e_to_g_super_enhancers = ps.sp_e_to_g_super_enhancers
        File? sp_g_to_e_super_enhancers = ps.sp_g_to_e_super_enhancers
        File? flankbedfile = if paired_sample_m then ps.flankbedfile else ss.flankbedfile
        File? ame_tsv = if paired_sample_m then ps.ame_tsv else ss.ame_tsv
        File? ame_html = if paired_sample_m then ps.ame_html else ss.ame_html
        File? ame_seq = if paired_sample_m then ps.ame_seq else ss.ame_seq
        File? meme = if paired_sample_m then ps.meme else ss.meme
        File? meme_summary = if paired_sample_m then ps.meme_summary else ss.meme_summary
        File? summit_ame_tsv = if paired_sample_m then ps.summit_ame_tsv else ss.summit_ame_tsv
        File? summit_ame_html = if paired_sample_m then ps.summit_ame_html else ss.summit_ame_html
        File? summit_ame_seq = if paired_sample_m then ps.summit_ame_seq else ss.summit_ame_seq
        File? summit_meme = if paired_sample_m then ps.summit_meme else ss.summit_meme
        File? summit_meme_summary = if paired_sample_m then ps.summit_meme_summary else ss.summit_meme_summary
        File? sp_flankbedfile = ps.summit_meme_summary
        File? sp_ame_tsv = ps.sp_ame_tsv
        File? sp_ame_html = ps.sp_ame_html
        File? sp_ame_seq = ps.sp_ame_seq
        File? sp_meme = ps.sp_meme
        File? sp_meme_summary = ps.sp_meme_summary
        File? sp_summit_ame_tsv = ps.sp_summit_ame_tsv
        File? sp_summit_ame_html = ps.sp_summit_ame_html
        File? sp_summit_ame_seq = ps.sp_summit_ame_seq
        File? sp_summit_meme = ps.sp_summit_meme
        File? sp_summit_meme_summary = ps.sp_summit_meme_summary
        File? s_matrices = if paired_sample_m then ps.s_matrices else ss.s_matrices
        File? densityplot = if paired_sample_m then ps.densityplot else ss.densityplot
        File? pdf_gene = if paired_sample_m then ps.pdf_gene else ss.pdf_gene
        File? pdf_h_gene = if paired_sample_m then ps.pdf_h_gene else ss.pdf_h_gene
        File? png_h_gene = if paired_sample_m then ps.png_h_gene else ss.png_h_gene
        File? jpg_h_gene = if paired_sample_m then ps.jpg_h_gene else ss.jpg_h_gene
        File? pdf_promoters = if paired_sample_m then ps.pdf_promoters else ss.pdf_promoters
        File? pdf_h_promoters = if paired_sample_m then ps.pdf_h_promoters else ss.pdf_h_promoters
        File? png_h_promoters = if paired_sample_m then ps.png_h_promoters else ss.png_h_promoters
        File? jpg_h_promoters = if paired_sample_m then ps.jpg_h_promoters else ss.jpg_h_promoters
        File? sp_s_matrices = ps.sp_s_matrices
        File? sp_densityplot = ps.sp_densityplot
        File? sp_pdf_gene = ps.sp_pdf_gene
        File? sp_pdf_h_gene = ps.sp_pdf_h_gene
        File? sp_png_h_gene = ps.sp_png_h_gene
        File? sp_jpg_h_gene = ps.sp_jpg_h_gene
        File? sp_pdf_promoters = ps.sp_pdf_promoters
        File? sp_pdf_h_promoters = ps.sp_pdf_h_promoters
        File? sp_png_h_promoters = ps.sp_png_h_promoters
        File? sp_jpg_h_promoters = ps.sp_jpg_h_promoters
        File? peak_promoters = if paired_sample_m then ps.peak_promoters else ss.peak_promoters
        File? peak_genebody = if paired_sample_m then ps.peak_genebody else ss.peak_genebody
        File? peak_window = if paired_sample_m then ps.peak_window else ss.peak_window
        File? peak_closest = if paired_sample_m then ps.peak_closest else ss.peak_closest
        File? peak_comparison = if paired_sample_m then ps.peak_comparison else ss.peak_comparison
        File? gene_comparison = if paired_sample_m then ps.gene_comparison else ss.gene_comparison
        File? pdf_comparison = if paired_sample_m then ps.pdf_comparison else ss.pdf_comparison
        File? all_peak_promoters = if paired_sample_m then ps.all_peak_promoters else ss.all_peak_promoters
        File? all_peak_genebody = if paired_sample_m then ps.all_peak_genebody else ss.all_peak_genebody
        File? all_peak_window = if paired_sample_m then ps.all_peak_window else ss.all_peak_window
        File? all_peak_closest = if paired_sample_m then ps.all_peak_closest else ss.all_peak_closest
        File? all_peak_comparison = if paired_sample_m then ps.all_peak_comparison else ss.all_peak_comparison
        File? all_gene_comparison = if paired_sample_m then ps.all_gene_comparison else ss.all_gene_comparison
        File? all_pdf_comparison = if paired_sample_m then ps.all_pdf_comparison else ss.all_pdf_comparison
        File? nomodel_peak_promoters = if paired_sample_m then ps.nomodel_peak_promoters else ss.nomodel_peak_promoters
        File? nomodel_peak_genebody = if paired_sample_m then ps.nomodel_peak_genebody else ss.nomodel_peak_genebody
        File? nomodel_peak_window = if paired_sample_m then ps.nomodel_peak_window else ss.nomodel_peak_window
        File? nomodel_peak_closest = if paired_sample_m then ps.nomodel_peak_closest else ss.nomodel_peak_closest
        File? nomodel_peak_comparison = if paired_sample_m then ps.nomodel_peak_comparison else ss.nomodel_peak_comparison
        File? nomodel_gene_comparison = if paired_sample_m then ps.nomodel_gene_comparison else ss.nomodel_gene_comparison
        File? nomodel_pdf_comparison = if paired_sample_m then ps.nomodel_pdf_comparison else ss.nomodel_pdf_comparison
        File? sicer_peak_promoters = if paired_sample_m then ps.sicer_peak_promoters else ss.sicer_peak_promoters
        File? sicer_peak_genebody = if paired_sample_m then ps.sicer_peak_genebody else ss.sicer_peak_genebody
        File? sicer_peak_window = if paired_sample_m then ps.sicer_peak_window else ss.sicer_peak_window
        File? sicer_peak_closest = if paired_sample_m then ps.sicer_peak_closest else ss.sicer_peak_closest
        File? sicer_peak_comparison = if paired_sample_m then ps.sicer_peak_comparison else ss.sicer_peak_comparison
        File? sicer_gene_comparison = if paired_sample_m then ps.sicer_gene_comparison else ss.sicer_gene_comparison
        File? sicer_pdf_comparison = if paired_sample_m then ps.sicer_pdf_comparison else ss.sicer_pdf_comparison
        File? sp_peak_promoters = ps.sp_peak_promoters
        File? sp_peak_genebody =  ps.sp_peak_genebody
        File? sp_peak_window = ps.sp_peak_window
        File? sp_peak_closest = ps.sp_peak_closest
        File? sp_peak_comparison = ps.sp_peak_comparison
        File? sp_gene_comparison = ps.sp_gene_comparison
        File? sp_pdf_comparison = ps.sp_pdf_comparison
        File? sp_all_peak_promoters = ps.sp_all_peak_promoters
        File? sp_all_peak_genebody =  ps.sp_all_peak_genebody
        File? sp_all_peak_window = ps.sp_all_peak_window
        File? sp_all_peak_closest = ps.sp_all_peak_closest
        File? sp_all_peak_comparison = ps.sp_all_peak_comparison
        File? sp_all_gene_comparison = ps.sp_all_gene_comparison
        File? sp_all_pdf_comparison = ps.sp_all_pdf_comparison
        File? sp_nomodel_peak_promoters = ps.sp_nomodel_peak_promoters
        File? sp_nomodel_peak_genebody =  ps.sp_nomodel_peak_genebody
        File? sp_nomodel_peak_window = ps.sp_nomodel_peak_window
        File? sp_nomodel_peak_closest = ps.sp_nomodel_peak_closest
        File? sp_nomodel_peak_comparison = ps.sp_nomodel_peak_comparison
        File? sp_nomodel_gene_comparison = ps.sp_nomodel_gene_comparison
        File? sp_nomodel_pdf_comparison = ps.sp_nomodel_pdf_comparison
        File? sp_sicer_peak_promoters = ps.sp_sicer_peak_promoters
        File? sp_sicer_peak_genebody =  ps.sp_sicer_peak_genebody
        File? sp_sicer_peak_window = ps.sp_sicer_peak_window
        File? sp_sicer_peak_closest = ps.sp_sicer_peak_closest
        File? sp_sicer_peak_comparison = ps.sp_sicer_peak_comparison
        File? sp_sicer_gene_comparison = ps.sp_sicer_gene_comparison
        File? sp_sicer_pdf_comparison = ps.sp_sicer_pdf_comparison
        File? bigwig = if paired_sample_m then ps.bigwig else ss.bigwig
        File? norm_wig = if paired_sample_m then ps.norm_wig else ss.norm_wig
        File? tdffile = if paired_sample_m then ps.tdffile else ss.tdffile
        File? n_bigwig = if paired_sample_m then ps.n_bigwig else ss.n_bigwig
        File? n_norm_wig = if paired_sample_m then ps.n_norm_wig else ss.n_norm_wig
        File? n_tdffile = if paired_sample_m then ps.n_tdffile else ss.n_tdffile
        File? a_bigwig = if paired_sample_m then ps.a_bigwig else ss.a_bigwig
        File? a_norm_wig = if paired_sample_m then ps.a_norm_wig else ss.a_norm_wig
        File? a_tdffile = if paired_sample_m then ps.a_tdffile else ss.a_tdffile
        File? s_bigwig = if paired_sample_m then ps.s_bigwig else ss.s_bigwig
        File? s_norm_wig = if paired_sample_m then ps.s_norm_wig else ss.s_norm_wig
        File? s_tdffile = if paired_sample_m then ps.s_tdffile else ss.s_tdffile
        File? sp_bigwig = ps.sp_bigwig
        File? sp_norm_wig = ps.sp_norm_wig
        File? sp_tdffile = ps.sp_tdffile
        File? sp_n_bigwig = ps.sp_n_bigwig
        File? sp_n_norm_wig = ps.sp_n_norm_wig
        File? sp_n_tdffile = ps.sp_n_tdffile
        File? sp_a_bigwig = ps.sp_a_bigwig
        File? sp_a_norm_wig = ps.sp_a_norm_wig
        File? sp_a_tdffile = ps.sp_a_tdffile
        File? sp_s_bigwig = ps.sp_s_bigwig
        File? sp_s_norm_wig = ps.sp_s_norm_wig
        File? sp_s_tdffile = ps.sp_s_tdffile
        File? sf_bigwig = ps.sf_bigwig
        File? sf_tdffile = ps.sf_tdffile
        File? sf_wigfile = ps.sf_wigfile
        Array[File?]? s_qc_statsfile = if paired_sample_m then ps.s_qc_statsfile else ss.s_qc_statsfile
        Array[File?]? s_qc_htmlfile = if paired_sample_m then ps.s_qc_htmlfile else ss.s_qc_htmlfile
        Array[File?]? s_qc_textfile = if paired_sample_m then ps.s_qc_textfile else ss.s_qc_textfile
        File? s_qc_mergehtml = if paired_sample_m then ps.s_qc_mergehtml else ss.s_qc_mergehtml
        Array[File?]? sp_qc_statsfile = ps.sp_qc_statsfile
        Array[File?]? sp_qc_htmlfile = ps.sp_qc_htmlfile
        Array[File?]? sp_qc_textfile = ps.sp_qc_textfile
        File? s_uno_statsfile = if paired_sample_m then ps.s_uno_statsfile else ss.s_uno_statsfile
        File? s_uno_htmlfile = if paired_sample_m then ps.s_uno_htmlfile else ss.s_uno_htmlfile
        File? s_uno_textfile = if paired_sample_m then ps.s_uno_textfile else ss.s_uno_textfile
        File? statsfile = if paired_sample_m then ps.statsfile else ss.statsfile
        File? htmlfile = if paired_sample_m then ps.htmlfile else ss.htmlfile
        File? textfile = if paired_sample_m then ps.textfile else ss.textfile
        File? sp_statsfile = ps.sp_statsfile
        File? sp_htmlfile = ps.sp_htmlfile
        File? sp_textfile = ps.sp_textfile
        File? s_summaryhtml = ps.s_summaryhtml
        File? s_summarytxt = ps.s_summarytxt
        File? s_fragsize = ps.s_fragsize
        File? summaryhtml = if paired_sample_m then ps.summaryhtml else ss.summaryhtml
        File? summarytxt = if paired_sample_m then ps.summarytxt else ss.summarytxt
        File? sp_summaryhtml = ps.sp_summaryhtml
        File? sp_summarytxt = ps.sp_summarytxt
    }
}
