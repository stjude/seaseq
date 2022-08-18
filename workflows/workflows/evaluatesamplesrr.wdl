version 1.0
import "../../seaseq-case.wdl" as ss
import "../../peaseq-case.wdl" as ps
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
        String? results_name
        Boolean run_motifs=true
        Int? insertsize = 600
        String? strandedness = "fr"
    }
    
    Boolean paired=true
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
                blacklist=blacklist,
                gtf=gtf,
                bowtie_index=bowtie_index,
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
                blacklist=blacklist,
                gtf=gtf,
                bowtie_index=bowtie_index,
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
        Array[File?]? indv_s_htmlfile = select_first([ss.indv_s_htmlfile, ps.indv_s_htmlfile])
        Array[File?]? indv_s_zipfile = select_first([ss.indv_s_zipfile, ps.indv_s_zipfile])
        Array[File?]? indv_s_bam_htmlfile = select_first([ss.indv_s_bam_htmlfile, ps.indv_s_bam_htmlfile])
        Array[File?]? indv_s_bam_zipfile = select_first([ss.indv_s_bam_zipfile, ps.indv_s_bam_zipfile])
        File? s_mergebam_htmlfile = select_first([ss.s_mergebam_htmlfile, ps.s_mergebam_htmlfile])
        File? s_mergebam_zipfile = select_first([ss.s_mergebam_zipfile, ps.s_mergebam_zipfile])
        Array[File?]? indv_sp_bam_htmlfile = ps.indv_sp_bam_htmlfile
        Array[File?]? indv_sp_bam_zipfile = ps.indv_sp_bam_zipfile
        File? sp_mergebam_htmlfile = ps.sp_mergebam_htmlfile
        File? sp_mergebam_zipfile = ps.sp_mergebam_zipfile
        File? uno_s_htmlfile = ss.uno_s_htmlfile
        File? uno_s_zipfile = ss.uno_s_zipfile
        File? uno_s_bam_htmlfile = select_first([ss.uno_s_bam_htmlfile, ps.uno_s_bam_htmlfile])
        File? uno_s_bam_zipfile = select_first([ss.uno_s_bam_zipfile, ps.uno_s_bam_zipfile])
        Array[File?]? s_metrics_out = select_first([ss.s_metrics_out, ps.s_metrics_out])
        File? uno_s_metrics_out = ss.uno_s_metrics_out
        Array[File?]? indv_s_sortedbam = select_first([ss.indv_s_sortedbam, ps.indv_s_sortedbam])
        Array[File?]? indv_s_indexbam = select_first([ss.indv_s_indexbam, ps.indv_s_indexbam])
        Array[File?]? indv_s_bkbam = select_first([ss.indv_s_bkbam, ps.indv_s_bkbam])
        Array[File?]? indv_s_bkindexbam = select_first([ss.indv_s_bkindexbam, ps.indv_s_bkindexbam])
        Array[File?]? indv_s_rmbam = select_first([ss.indv_s_rmbam, ps.indv_s_rmbam])
        Array[File?]? indv_s_rmindexbam = select_first([ss.indv_s_rmindexbam, ps.indv_s_rmindexbam])
        Array[File?]? indv_sp_sortedbam = ps.indv_sp_sortedbam
        Array[File?]? indv_sp_indexbam = ps.indv_sp_indexbam
        Array[File?]? indv_sp_bkbam = ps.indv_sp_bkbam
        Array[File?]? indv_sp_bkindexbam = ps.indv_sp_bkindexbam
        Array[File?]? indv_sp_rmbam = ps.indv_sp_rmbam
        Array[File?]? indv_sp_rmindexbam = ps.indv_sp_rmindexbam
        File? uno_s_sortedbam = select_first([ss.uno_s_sortedbam, ps.uno_s_sortedbam])
        File? uno_s_indexstatsbam = select_first([ss.uno_s_indexstatsbam, ps.uno_s_indexstatsbam])
        File? uno_s_bkbam = select_first([ss.uno_s_bkbam, ps.uno_s_bkbam])
        File? uno_s_bkindexbam = select_first([ss.uno_s_bkindexbam, ps.uno_s_bkindexbam])
        File? uno_s_rmbam = select_first([ss.uno_s_rmbam, ps.uno_s_rmbam])
        File? uno_s_rmindexbam = select_first([ss.uno_s_rmindexbam, ps.uno_s_rmindexbam])
        File? s_mergebamfile = select_first([ss.s_mergebamfile, ps.s_mergebamfile])
        File? s_mergebamindex = select_first([ss.s_mergebamindex, ps.s_mergebamindex])
        File? s_bkbam = select_first([ss.s_bkbam, ps.s_bkbam])
        File? s_bkindexbam = select_first([ss.s_bkindexbam, ps.s_bkindexbam])
        File? s_rmbam = select_first([ss.s_rmbam, ps.s_rmbam])
        File? s_rmindexbam = select_first([ss.s_rmindexbam, ps.s_rmindexbam])
        File? sp_mergebamfile = ps.sp_mergebamfile
        File? sp_mergebamindex = ps.sp_mergebamindex
        File? sp_bkbam = ps.sp_bkbam
        File? sp_bkindexbam = ps.sp_bkindexbam
        File? sp_rmbam = ps.sp_rmbam
        File? sp_rmindexbam = ps.sp_rmindexbam
        File? s_fragments_bam = ps.s_fragments_bam
        File? s_fragments_indexbam = ps.s_fragments_indexbam
        File? peakbedfile = select_first([ss.peakbedfile, ps.peakbedfile])
        File? peakxlsfile = select_first([ss.peakxlsfile, ps.peakxlsfile])
        File? negativexlsfile = select_first([ss.negativexlsfile, ps.negativexlsfile])
        File? summitsfile = select_first([ss.summitsfile, ps.summitsfile])
        File? wigfile = select_first([ss.wigfile, ps.wigfile])
        File? all_peakbedfile = select_first([ss.all_peakbedfile, ps.all_peakbedfile])
        File? all_peakxlsfile = select_first([ss.all_peakxlsfile, ps.all_peakxlsfile])
        File? all_negativexlsfile = select_first([ss.all_negativexlsfile, ps.all_negativexlsfile])
        File? all_summitsfile = select_first([ss.all_summitsfile, ps.all_summitsfile])
        File? all_wigfile = select_first([ss.all_wigfile, ps.all_wigfile])
        File? nm_peakbedfile = select_first([ss.nm_peakbedfile, ps.nm_peakbedfile])
        File? nm_peakxlsfile = select_first([ss.nm_peakxlsfile, ps.nm_peakxlsfile])
        File? nm_negativexlsfile = select_first([ss.nm_negativexlsfile, ps.nm_negativexlsfile])
        File? nm_summitsfile = select_first([ss.nm_summitsfile, ps.nm_summitsfile])
        File? nm_wigfile = select_first([ss.nm_wigfile, ps.nm_wigfile])
        File? readme_peaks = select_first([ss.readme_peaks, ps.readme_peaks])
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
        File? scoreisland = select_first([ss.scoreisland, ps.scoreisland])
        File? sicer_wigfile = select_first([ss.sicer_wigfile, ps.sicer_wigfile])
        File? sp_scoreisland = ps.sp_scoreisland
        File? sp_sicer_wigfile = ps.sp_sicer_wigfile
        File? pngfile = select_first([ss.pngfile, ps.pngfile])
        File? mapped_union = select_first([ss.mapped_union, ps.mapped_union])
        File? mapped_stitch = select_first([ss.mapped_stitch, ps.mapped_stitch])
        File? enhancers = select_first([ss.enhancers, ps.enhancers])
        File? super_enhancers = select_first([ss.super_enhancers, ps.super_enhancers])
        File? gff_file = select_first([ss.gff_file, ps.gff_file])
        File? gff_union = select_first([ss.gff_union, ps.gff_union])
        File? union_enhancers = select_first([ss.union_enhancers, ps.union_enhancers])
        File? stitch_enhancers = select_first([ss.stitch_enhancers, ps.stitch_enhancers])
        File? e_to_g_enhancers = select_first([ss.e_to_g_enhancers, ps.e_to_g_enhancers])
        File? g_to_e_enhancers = select_first([ss.g_to_e_enhancers, ps.g_to_e_enhancers])
        File? e_to_g_super_enhancers = select_first([ss.e_to_g_super_enhancers, ps.e_to_g_super_enhancers])
        File? g_to_e_super_enhancers = select_first([ss.g_to_e_super_enhancers, ps.g_to_e_super_enhancers])
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
        File? flankbedfile = select_first([ss.flankbedfile, ps.flankbedfile])
        File? ame_tsv = select_first([ss.ame_tsv, ps.ame_tsv])
        File? ame_html = select_first([ss.ame_html, ps.ame_html])
        File? ame_seq = select_first([ss.ame_seq, ps.ame_seq])
        File? meme = select_first([ss.meme, ps.meme])
        File? meme_summary = select_first([ss.meme_summary, ps.meme_summary])
        File? summit_ame_tsv = select_first([ss.summit_ame_tsv, ps.summit_ame_tsv])
        File? summit_ame_html = select_first([ss.summit_ame_html, ps.summit_ame_html])
        File? summit_ame_seq = select_first([ss.summit_ame_seq, ps.summit_ame_seq])
        File? summit_meme = select_first([ss.summit_meme, ps.summit_meme])
        File? summit_meme_summary = select_first([ss.summit_meme_summary, ps.summit_meme_summary])
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
        File? s_matrices = select_first([ss.s_matrices, ps.s_matrices])
        File? densityplot = select_first([ss.densityplot, ps.densityplot])
        File? pdf_gene = select_first([ss.pdf_gene, ps.pdf_gene])
        File? pdf_h_gene = select_first([ss.pdf_h_gene, ps.pdf_h_gene])
        File? png_h_gene = select_first([ss.png_h_gene, ps.png_h_gene])
        File? jpg_h_gene = select_first([ss.jpg_h_gene, ps.jpg_h_gene])
        File? pdf_promoters = select_first([ss.pdf_promoters, ps.pdf_promoters])
        File? pdf_h_promoters = select_first([ss.pdf_h_promoters, ps.pdf_h_promoters])
        File? png_h_promoters = select_first([ss.png_h_promoters, ps.png_h_promoters])
        File? jpg_h_promoters = select_first([ss.jpg_h_promoters, ps.jpg_h_promoters])
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
        File? peak_promoters = select_first([ss.peak_promoters, ps.peak_promoters])
        File? peak_genebody = select_first([ss.peak_genebody, ps.peak_genebody])
        File? peak_window = select_first([ss.peak_window, ps.peak_window])
        File? peak_closest = select_first([ss.peak_closest, ps.peak_closest])
        File? peak_comparison = select_first([ss.peak_comparison, ps.peak_comparison])
        File? gene_comparison = select_first([ss.gene_comparison, ps.gene_comparison])
        File? pdf_comparison = select_first([ss.pdf_comparison, ps.pdf_comparison])
        File? all_peak_promoters = select_first([ss.all_peak_promoters, ps.all_peak_promoters])
        File? all_peak_genebody = select_first([ss.all_peak_genebody, ps.all_peak_genebody])
        File? all_peak_window = select_first([ss.all_peak_window, ps.all_peak_window])
        File? all_peak_closest = select_first([ss.all_peak_closest, ps.all_peak_closest])
        File? all_peak_comparison = select_first([ss.all_peak_comparison, ps.all_peak_comparison])
        File? all_gene_comparison = select_first([ss.all_gene_comparison, ps.all_gene_comparison])
        File? all_pdf_comparison = select_first([ss.all_pdf_comparison, ps.all_pdf_comparison])
        File? nomodel_peak_promoters = select_first([ss.nomodel_peak_promoters, ps.nomodel_peak_promoters])
        File? nomodel_peak_genebody = select_first([ss.nomodel_peak_genebody, ps.nomodel_peak_genebody])
        File? nomodel_peak_window = select_first([ss.nomodel_peak_window, ps.nomodel_peak_window])
        File? nomodel_peak_closest = select_first([ ss.nomodel_peak_closest, ps.nomodel_peak_closest])
        File? nomodel_peak_comparison = select_first([ss.nomodel_peak_comparison, ps.nomodel_peak_comparison])
        File? nomodel_gene_comparison = select_first([ss.nomodel_gene_comparison, ps.nomodel_gene_comparison])
        File? nomodel_pdf_comparison = select_first([ss.nomodel_pdf_comparison, ps.nomodel_pdf_comparison])
        File? sicer_peak_promoters = select_first([ss.sicer_peak_promoters, ps.sicer_peak_promoters])
        File? sicer_peak_genebody = select_first([ss.sicer_peak_genebody, ps.sicer_peak_genebody])
        File? sicer_peak_window = select_first([ss.sicer_peak_window, ps.sicer_peak_window])
        File? sicer_peak_closest = select_first([ss.sicer_peak_closest, ps.sicer_peak_closest])
        File? sicer_peak_comparison = select_first([ss.sicer_peak_comparison, ps.sicer_peak_comparison])
        File? sicer_gene_comparison = select_first([ss.sicer_gene_comparison, ps.sicer_gene_comparison])
        File? sicer_pdf_comparison = select_first([ss.sicer_pdf_comparison, ps.sicer_pdf_comparison])
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
        File? bigwig = select_first([ss.bigwig, ps.bigwig])
        File? norm_wig = select_first([ss.norm_wig, ps.norm_wig])
        File? tdffile = select_first([ss.tdffile, ps.tdffile])
        File? n_bigwig = select_first([ss.n_bigwig, ps.n_bigwig])
        File? n_norm_wig = select_first([ss.n_norm_wig, ps.n_norm_wig])
        File? n_tdffile = select_first([ss.n_tdffile, ps.n_tdffile])
        File? a_bigwig = select_first([ss.a_bigwig, ps.a_bigwig])
        File? a_norm_wig = select_first([ss.a_norm_wig, ps.a_norm_wig])
        File? a_tdffile = select_first([ss.a_tdffile, ps.a_tdffile])
        File? s_bigwig = select_first([ss.s_bigwig, ps.s_bigwig])
        File? s_norm_wig = select_first([ss.s_norm_wig, ps.s_norm_wig])
        File? s_tdffile = select_first([ss.s_tdffile, ps.s_tdffile])
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
        Array[File?]? s_qc_statsfile = select_first([ss.s_qc_statsfile, ps.s_qc_statsfile])
        Array[File?]? s_qc_htmlfile = select_first([ss.s_qc_htmlfile, ps.s_qc_htmlfile])
        Array[File?]? s_qc_textfile = select_first([ss.s_qc_textfile, ps.s_qc_textfile])
        File? s_qc_mergehtml = select_first([ss.s_qc_mergehtml, ps.s_qc_mergehtml])
        Array[File?]? sp_qc_statsfile = ps.sp_qc_statsfile
        Array[File?]? sp_qc_htmlfile = ps.sp_qc_htmlfile
        Array[File?]? sp_qc_textfile = ps.sp_qc_textfile
        File? s_uno_statsfile = select_first([ss.s_uno_statsfile, ps.s_uno_statsfile])
        File? s_uno_htmlfile = select_first([ss.s_uno_htmlfile, ps.s_uno_htmlfile])
        File? s_uno_textfile = select_first([ss.s_uno_textfile, ps.s_uno_textfile])
        File? statsfile = select_first([ss.statsfile, ps.statsfile])
        File? htmlfile = select_first([ss.htmlfile, ps.htmlfile])
        File? textfile = select_first([ss.textfile, ps.textfile])
        File? sp_statsfile = ps.sp_statsfile
        File? sp_htmlfile = ps.sp_htmlfile
        File? sp_textfile = ps.sp_textfile
        File? s_summaryhtml = ps.s_summaryhtml
        File? s_summarytxt = ps.s_summarytxt
        File? s_fragsize = ps.s_fragsize
        File? summaryhtml = select_first([ss.summaryhtml, ps.summaryhtml])
        File? summarytxt = select_first([ss.summarytxt, ps.summarytxt])
        File? sp_summaryhtml = ps.sp_summaryhtml
        File? sp_summarytxt = ps.sp_summarytxt
    }
}
