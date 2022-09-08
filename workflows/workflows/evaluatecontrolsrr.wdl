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
        Boolean paired=true
    }
    
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
        Array[File?]? indv_s_htmlfile = if (paired_sample_m && paired_control_m) then pc.indv_s_htmlfile else sc.indv_s_htmlfile
        Array[File?]? indv_s_zipfile = if (paired_sample_m && paired_control_m) then pc.indv_s_zipfile else sc.indv_s_zipfile
        Array[File?]? indv_s_bam_htmlfile = if (paired_sample_m && paired_control_m) then pc.indv_s_bam_htmlfile else sc.indv_s_bam_htmlfile
        Array[File?]? indv_s_bam_zipfile = if (paired_sample_m && paired_control_m) then pc.indv_s_bam_zipfile else sc.indv_s_bam_zipfile
        Array[File?]? indv_c_htmlfile = if (paired_sample_m && paired_control_m) then pc.indv_c_htmlfile else sc.indv_c_htmlfile
        Array[File?]? indv_c_zipfile = if (paired_sample_m && paired_control_m) then pc.indv_c_zipfile else sc.indv_c_zipfile
        Array[File?]? indv_c_bam_htmlfile = if (paired_sample_m && paired_control_m) then pc.indv_c_bam_htmlfile else sc.indv_c_bam_htmlfile
        Array[File?]? indv_c_bam_zipfile = if (paired_sample_m && paired_control_m) then pc.indv_c_bam_zipfile else sc.indv_c_bam_zipfile
        File? s_mergebam_htmlfile = if (paired_sample_m && paired_control_m) then pc.s_mergebam_htmlfile else sc.s_mergebam_htmlfile
        File? s_mergebam_zipfile = if (paired_sample_m && paired_control_m) then pc.s_mergebam_zipfile else sc.s_mergebam_zipfile
        File? c_mergebam_htmlfile = if (paired_sample_m && paired_control_m) then pc.c_mergebam_htmlfile else sc.c_mergebam_htmlfile
        File? c_mergebam_zipfile = if (paired_sample_m && paired_control_m) then pc.c_mergebam_zipfile else sc.c_mergebam_zipfile
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
        File? uno_s_bam_htmlfile = if (paired_sample_m && paired_control_m) then pc.uno_s_bam_htmlfile else sc.uno_s_bam_htmlfile
        File? uno_s_bam_zipfile = if (paired_sample_m && paired_control_m) then pc.uno_s_bam_zipfile else sc.uno_s_bam_zipfile
        File? uno_c_htmlfile = sc.uno_c_htmlfile
        File? uno_c_zipfile = sc.uno_c_zipfile
        File? uno_c_bam_htmlfile = if (paired_sample_m && paired_control_m) then pc.uno_c_bam_htmlfile else sc.uno_c_bam_htmlfile
        File? uno_c_bam_zipfile = if (paired_sample_m && paired_control_m) then pc.uno_c_bam_zipfile else sc.uno_c_bam_zipfile
        Array[File?]? s_metrics_out = if (paired_sample_m && paired_control_m) then pc.s_metrics_out else sc.s_metrics_out
        File? uno_s_metrics_out = sc.uno_s_metrics_out
        Array[File?]? c_metrics_out = if (paired_sample_m && paired_control_m) then pc.c_metrics_out else sc.c_metrics_out
        File? uno_c_metrics_out = sc.uno_c_metrics_out
        Array[File?]? indv_s_sortedbam = if (paired_sample_m && paired_control_m) then pc.indv_s_sortedbam else sc.indv_s_sortedbam
        Array[File?]? indv_s_indexbam = if (paired_sample_m && paired_control_m) then pc.indv_s_indexbam else sc.indv_s_indexbam
        Array[File?]? indv_s_bkbam = if (paired_sample_m && paired_control_m) then pc.indv_s_bkbam else sc.indv_s_bkbam
        Array[File?]? indv_s_bkindexbam = if (paired_sample_m && paired_control_m) then pc.indv_s_bkindexbam else sc.indv_s_bkindexbam
        Array[File?]? indv_s_rmbam = if (paired_sample_m && paired_control_m) then pc.indv_s_rmbam else sc.indv_s_rmbam
        Array[File?]? indv_s_rmindexbam = if (paired_sample_m && paired_control_m) then pc.indv_s_rmindexbam else sc.indv_s_rmindexbam
        Array[File?]? indv_c_sortedbam = if (paired_sample_m && paired_control_m) then pc.indv_c_sortedbam else sc.indv_c_sortedbam
        Array[File?]? indv_c_indexbam = if (paired_sample_m && paired_control_m) then pc.indv_c_indexbam else sc.indv_c_indexbam
        Array[File?]? indv_c_bkbam = if (paired_sample_m && paired_control_m) then pc.indv_c_bkbam else sc.indv_c_bkbam
        Array[File?]? indv_c_bkindexbam = if (paired_sample_m && paired_control_m) then pc.indv_c_bkindexbam else sc.indv_c_bkindexbam
        Array[File?]? indv_c_rmbam = if (paired_sample_m && paired_control_m) then pc.indv_c_rmbam else sc.indv_c_rmbam
        Array[File?]? indv_c_rmindexbam = if (paired_sample_m && paired_control_m) then pc.indv_c_rmindexbam else sc.indv_c_rmindexbam
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
        File? uno_s_sortedbam = if (paired_sample_m && paired_control_m) then pc.uno_s_sortedbam else sc.uno_s_sortedbam
        File? uno_s_indexstatsbam = if (paired_sample_m && paired_control_m) then pc.uno_s_indexstatsbam else sc.uno_s_indexstatsbam
        File? uno_s_bkbam = if (paired_sample_m && paired_control_m) then pc.uno_s_bkbam else sc.uno_s_bkbam
        File? uno_s_bkindexbam = if (paired_sample_m && paired_control_m) then pc.uno_s_bkindexbam else sc.uno_s_bkindexbam
        File? uno_s_rmbam = if (paired_sample_m && paired_control_m) then pc.uno_s_rmbam else sc.uno_s_rmbam
        File? uno_s_rmindexbam = if (paired_sample_m && paired_control_m) then pc.uno_s_rmindexbam else sc.uno_s_rmindexbam
        File? uno_c_sortedbam = if (paired_sample_m && paired_control_m) then pc.uno_c_sortedbam else sc.uno_c_sortedbam
        File? uno_c_indexstatsbam = if (paired_sample_m && paired_control_m) then pc.uno_c_indexstatsbam else sc.uno_c_indexstatsbam
        File? uno_c_bkbam = if (paired_sample_m && paired_control_m) then pc.uno_c_bkbam else sc.uno_c_bkbam
        File? uno_c_bkindexbam = if (paired_sample_m && paired_control_m) then pc.uno_c_bkindexbam else sc.uno_c_bkindexbam
        File? uno_c_rmbam = if (paired_sample_m && paired_control_m) then pc.uno_c_rmbam else sc.uno_c_rmbam
        File? uno_c_rmindexbam = if (paired_sample_m && paired_control_m) then pc.uno_c_rmindexbam else sc.uno_c_rmindexbam
        File? s_mergebamfile = if (paired_sample_m && paired_control_m) then pc.s_mergebamfile else sc.s_mergebamfile
        File? s_mergebamindex = if (paired_sample_m && paired_control_m) then pc.s_mergebamindex else sc.s_mergebamindex
        File? s_bkbam = if (paired_sample_m && paired_control_m) then pc.s_bkbam else sc.s_bkbam
        File? s_bkindexbam = if (paired_sample_m && paired_control_m) then pc.s_bkindexbam else sc.s_bkindexbam
        File? s_rmbam = if (paired_sample_m && paired_control_m) then pc.s_rmbam else sc.s_rmbam
        File? s_rmindexbam = if (paired_sample_m && paired_control_m) then pc.s_rmindexbam else sc.s_rmindexbam
        File? c_mergebamfile = if (paired_sample_m && paired_control_m) then pc.c_mergebamfile else sc.c_mergebamfile
        File? c_mergebamindex = if (paired_sample_m && paired_control_m) then pc.c_mergebamindex else sc.c_mergebamindex
        File? c_bkbam = if (paired_sample_m && paired_control_m) then pc.c_bkbam else sc.c_bkbam
        File? c_bkindexbam = if (paired_sample_m && paired_control_m) then pc.c_bkindexbam else sc.c_bkindexbam
        File? c_rmbam = if (paired_sample_m && paired_control_m) then pc.c_rmbam else sc.c_rmbam
        File? c_rmindexbam = if (paired_sample_m && paired_control_m) then pc.c_rmindexbam else sc.c_rmindexbam
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
        File? peakbedfile = if (paired_sample_m && paired_control_m) then pc.peakbedfile else sc.peakbedfile
        File? peakxlsfile = if (paired_sample_m && paired_control_m) then pc.peakxlsfile else sc.peakxlsfile
        File? negativexlsfile = if (paired_sample_m && paired_control_m) then pc.negativexlsfile else sc.negativexlsfile
        File? summitsfile = if (paired_sample_m && paired_control_m) then pc.summitsfile else sc.summitsfile
        File? wigfile = if (paired_sample_m && paired_control_m) then pc.wigfile else sc.wigfile
        File? ctrlwigfile = if (paired_sample_m && paired_control_m) then pc.ctrlwigfile else sc.ctrlwigfile
        File? all_peakbedfile = if (paired_sample_m && paired_control_m) then pc.all_peakbedfile else sc.all_peakbedfile
        File? all_peakxlsfile = if (paired_sample_m && paired_control_m) then pc.all_peakxlsfile else sc.all_peakxlsfile
        File? all_negativexlsfile = if (paired_sample_m && paired_control_m) then pc.all_negativexlsfile else sc.all_negativexlsfile
        File? all_summitsfile = if (paired_sample_m && paired_control_m) then pc.all_summitsfile else sc.all_summitsfile
        File? all_wigfile = if (paired_sample_m && paired_control_m) then pc.all_wigfile else sc.all_wigfile
        File? all_ctrlwigfile = if (paired_sample_m && paired_control_m) then pc.all_ctrlwigfile else sc.all_ctrlwigfile
        File? nm_peakbedfile = if (paired_sample_m && paired_control_m) then pc.nm_peakbedfile else sc.nm_peakbedfile
        File? nm_peakxlsfile = if (paired_sample_m && paired_control_m) then pc.nm_peakxlsfile else sc.nm_peakxlsfile
        File? nm_negativexlsfile = if (paired_sample_m && paired_control_m) then pc.nm_negativexlsfile else sc.nm_negativexlsfile
        File? nm_summitsfile = if (paired_sample_m && paired_control_m) then pc.nm_summitsfile else sc.nm_summitsfile
        File? nm_wigfile = if (paired_sample_m && paired_control_m) then pc.nm_wigfile else sc.nm_wigfile
        File? nm_ctrlwigfile = if (paired_sample_m && paired_control_m) then pc.nm_ctrlwigfile else sc.nm_ctrlwigfile
        File? readme_peaks = if (paired_sample_m && paired_control_m) then pc.readme_peaks else sc.readme_peaks
        File? only_s_peakbedfile = if (paired_sample_m && paired_control_m) then pc.only_s_peakbedfile else sc.only_s_peakbedfile
        File? only_s_peakxlsfile = if (paired_sample_m && paired_control_m) then pc.only_s_peakxlsfile else sc.only_s_peakxlsfile
        File? only_s_summitsfile = if (paired_sample_m && paired_control_m) then pc.only_s_summitsfile else sc.only_s_summitsfile
        File? only_s_wigfile = if (paired_sample_m && paired_control_m) then pc.only_s_wigfile else sc.only_s_wigfile
        File? only_c_peakbedfile = if (paired_sample_m && paired_control_m) then pc.only_c_peakbedfile else sc.only_c_peakbedfile
        File? only_c_peakxlsfile = if (paired_sample_m && paired_control_m) then pc.only_c_peakxlsfile else sc.only_c_peakxlsfile
        File? only_c_summitsfile = if (paired_sample_m && paired_control_m) then pc.only_c_summitsfile else sc.only_c_summitsfile
        File? only_c_wigfile = if (paired_sample_m && paired_control_m) then pc.only_c_wigfile else sc.only_c_wigfile
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
        File? scoreisland = if (paired_sample_m && paired_control_m) then pc.scoreisland else sc.scoreisland
        File? sicer_wigfile = if (paired_sample_m && paired_control_m) then pc.sicer_wigfile else sc.sicer_wigfile
        File? sicer_summary = if (paired_sample_m && paired_control_m) then pc.sicer_summary else sc.sicer_summary
        File? sicer_fdrisland = if (paired_sample_m && paired_control_m) then pc.sicer_fdrisland else sc.sicer_fdrisland
        File? sp_scoreisland = pc.sp_scoreisland
        File? sp_sicer_wigfile = pc.sp_sicer_wigfile
        File? sp_sicer_summary = pc.sp_sicer_summary
        File? sp_sicer_fdrisland = pc.sp_sicer_fdrisland
        File? pngfile = if (paired_sample_m && paired_control_m) then pc.pngfile else sc.pngfile
        File? mapped_union = if (paired_sample_m && paired_control_m) then pc.mapped_union else sc.mapped_union
        File? mapped_stitch = if (paired_sample_m && paired_control_m) then pc.mapped_stitch else sc.mapped_stitch
        File? enhancers = if (paired_sample_m && paired_control_m) then pc.enhancers else sc.enhancers
        File? super_enhancers = if (paired_sample_m && paired_control_m) then pc.super_enhancers else sc.super_enhancers
        File? gff_file = if (paired_sample_m && paired_control_m) then pc.gff_file else sc.gff_file
        File? gff_union = if (paired_sample_m && paired_control_m) then pc.gff_union else sc.gff_union
        File? union_enhancers = if (paired_sample_m && paired_control_m) then pc.union_enhancers else sc.union_enhancers
        File? stitch_enhancers = if (paired_sample_m && paired_control_m) then pc.stitch_enhancers else sc.stitch_enhancers
        File? e_to_g_enhancers = if (paired_sample_m && paired_control_m) then pc.e_to_g_enhancers else sc.e_to_g_enhancers
        File? g_to_e_enhancers = if (paired_sample_m && paired_control_m) then pc.g_to_e_enhancers else sc.g_to_e_enhancers
        File? e_to_g_super_enhancers = if (paired_sample_m && paired_control_m) then pc.e_to_g_super_enhancers else sc.e_to_g_super_enhancers
        File? g_to_e_super_enhancers = if (paired_sample_m && paired_control_m) then pc.g_to_e_super_enhancers else sc.g_to_e_super_enhancers
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
        File? flankbedfile = if (paired_sample_m && paired_control_m) then pc.flankbedfile else sc.flankbedfile
        File? ame_tsv = if (paired_sample_m && paired_control_m) then pc.ame_tsv else sc.ame_tsv
        File? ame_html = if (paired_sample_m && paired_control_m) then pc.ame_html else sc.ame_html
        File? ame_seq = if (paired_sample_m && paired_control_m) then pc.ame_seq else sc.ame_seq
        File? meme = if (paired_sample_m && paired_control_m) then pc.meme else sc.meme
        File? meme_summary = if (paired_sample_m && paired_control_m) then pc.meme_summary else sc.meme_summary
        File? summit_ame_tsv = if (paired_sample_m && paired_control_m) then pc.summit_ame_tsv else sc.summit_ame_tsv
        File? summit_ame_html = if (paired_sample_m && paired_control_m) then pc.summit_ame_html else sc.summit_ame_html
        File? summit_ame_seq = if (paired_sample_m && paired_control_m) then pc.summit_ame_seq else sc.summit_ame_seq
        File? summit_meme = if (paired_sample_m && paired_control_m) then pc.summit_meme else sc.summit_meme
        File? summit_meme_summary = if (paired_sample_m && paired_control_m) then pc.summit_meme_summary else sc.summit_meme_summary
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
        File? s_matrices = if (paired_sample_m && paired_control_m) then pc.s_matrices else sc.s_matrices
        File? c_matrices = if (paired_sample_m && paired_control_m) then pc.c_matrices else sc.c_matrices
        File? densityplot = if (paired_sample_m && paired_control_m) then pc.densityplot else sc.densityplot
        File? pdf_gene = if (paired_sample_m && paired_control_m) then pc.pdf_gene else sc.pdf_gene
        File? pdf_h_gene = if (paired_sample_m && paired_control_m) then pc.pdf_h_gene else sc.pdf_h_gene
        File? png_h_gene = if (paired_sample_m && paired_control_m) then pc.png_h_gene else sc.png_h_gene
        File? jpg_h_gene = if (paired_sample_m && paired_control_m) then pc.jpg_h_gene else sc.jpg_h_gene
        File? pdf_promoters = if (paired_sample_m && paired_control_m) then pc.pdf_promoters else sc.pdf_promoters
        File? pdf_h_promoters = if (paired_sample_m && paired_control_m) then pc.pdf_h_promoters else sc.pdf_h_promoters
        File? png_h_promoters = if (paired_sample_m && paired_control_m) then pc.png_h_promoters else sc.png_h_promoters
        File? jpg_h_promoters = if (paired_sample_m && paired_control_m) then pc.jpg_h_promoters else sc.jpg_h_promoters
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
        File? peak_promoters = if (paired_sample_m && paired_control_m) then pc.peak_promoters else sc.peak_promoters
        File? peak_genebody = if (paired_sample_m && paired_control_m) then pc.peak_genebody else sc.peak_genebody
        File? peak_window = if (paired_sample_m && paired_control_m) then pc.peak_window else sc.peak_window
        File? peak_closest = if (paired_sample_m && paired_control_m) then pc.peak_closest else sc.peak_closest
        File? peak_comparison = if (paired_sample_m && paired_control_m) then pc.peak_comparison else sc.peak_comparison
        File? gene_comparison = if (paired_sample_m && paired_control_m) then pc.gene_comparison else sc.gene_comparison
        File? pdf_comparison = if (paired_sample_m && paired_control_m) then pc.pdf_comparison else sc.pdf_comparison
        File? all_peak_promoters = if (paired_sample_m && paired_control_m) then pc.all_peak_promoters else sc.all_peak_promoters
        File? all_peak_genebody = if (paired_sample_m && paired_control_m) then pc.all_peak_genebody else sc.all_peak_genebody
        File? all_peak_window = if (paired_sample_m && paired_control_m) then pc.all_peak_window else sc.all_peak_window
        File? all_peak_closest = if (paired_sample_m && paired_control_m) then pc.all_peak_closest else sc.all_peak_closest
        File? all_peak_comparison = if (paired_sample_m && paired_control_m) then pc.all_peak_comparison else sc.all_peak_comparison
        File? all_gene_comparison = if (paired_sample_m && paired_control_m) then pc.all_gene_comparison else sc.all_gene_comparison
        File? all_pdf_comparison = if (paired_sample_m && paired_control_m) then pc.all_pdf_comparison else sc.all_pdf_comparison
        File? nomodel_peak_promoters = if (paired_sample_m && paired_control_m) then pc.nomodel_peak_promoters else sc.nomodel_peak_promoters
        File? nomodel_peak_genebody = if (paired_sample_m && paired_control_m) then pc.nomodel_peak_genebody else sc.nomodel_peak_genebody
        File? nomodel_peak_window = if (paired_sample_m && paired_control_m) then pc.nomodel_peak_window else sc.nomodel_peak_window
        File? nomodel_peak_closest = if (paired_sample_m && paired_control_m) then pc.nomodel_peak_closest else sc.nomodel_peak_closest
        File? nomodel_peak_comparison = if (paired_sample_m && paired_control_m) then pc.nomodel_peak_comparison else sc.nomodel_peak_comparison
        File? nomodel_gene_comparison = if (paired_sample_m && paired_control_m) then pc.nomodel_gene_comparison else sc.nomodel_gene_comparison
        File? nomodel_pdf_comparison = if (paired_sample_m && paired_control_m) then pc.nomodel_pdf_comparison else sc.nomodel_pdf_comparison
        File? sicer_peak_promoters = if (paired_sample_m && paired_control_m) then pc.sicer_peak_promoters else sc.sicer_peak_promoters
        File? sicer_peak_genebody = if (paired_sample_m && paired_control_m) then pc.sicer_peak_genebody else sc.sicer_peak_genebody
        File? sicer_peak_window = if (paired_sample_m && paired_control_m) then pc.sicer_peak_window else sc.sicer_peak_window
        File? sicer_peak_closest = if (paired_sample_m && paired_control_m) then pc.sicer_peak_closest else sc.sicer_peak_closest
        File? sicer_peak_comparison = if (paired_sample_m && paired_control_m) then pc.sicer_peak_comparison else sc.sicer_peak_comparison
        File? sicer_gene_comparison = if (paired_sample_m && paired_control_m) then pc.sicer_gene_comparison else sc.sicer_gene_comparison
        File? sicer_pdf_comparison = if (paired_sample_m && paired_control_m) then pc.sicer_pdf_comparison else sc.sicer_pdf_comparison
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
        File? bigwig = if (paired_sample_m && paired_control_m) then pc.bigwig else sc.bigwig
        File? norm_wig = if (paired_sample_m && paired_control_m) then pc.norm_wig else sc.norm_wig
        File? tdffile = if (paired_sample_m && paired_control_m) then pc.tdffile else sc.tdffile
        File? n_bigwig = if (paired_sample_m && paired_control_m) then pc.n_bigwig else sc.n_bigwig
        File? n_norm_wig = if (paired_sample_m && paired_control_m) then pc.n_norm_wig else sc.n_norm_wig
        File? n_tdffile = if (paired_sample_m && paired_control_m) then pc.n_tdffile else sc.n_tdffile
        File? a_bigwig = if (paired_sample_m && paired_control_m) then pc.a_bigwig else sc.a_bigwig
        File? a_norm_wig = if (paired_sample_m && paired_control_m) then pc.a_norm_wig else sc.a_norm_wig
        File? a_tdffile = if (paired_sample_m && paired_control_m) then pc.a_tdffile else sc.a_tdffile
        File? c_bigwig = if (paired_sample_m && paired_control_m) then pc.c_bigwig else sc.c_bigwig
        File? c_norm_wig = if (paired_sample_m && paired_control_m) then pc.c_norm_wig else sc.c_norm_wig
        File? c_tdffile = if (paired_sample_m && paired_control_m) then pc.c_tdffile else sc.c_tdffile
        File? c_n_bigwig = if (paired_sample_m && paired_control_m) then pc.c_n_bigwig else sc.c_n_bigwig
        File? c_n_norm_wig = if (paired_sample_m && paired_control_m) then pc.c_n_norm_wig else sc.c_n_norm_wig
        File? c_n_tdffile = if (paired_sample_m && paired_control_m) then pc.c_n_tdffile else sc.c_n_tdffile
        File? c_a_bigwig = if (paired_sample_m && paired_control_m) then pc.c_a_bigwig else sc.c_a_bigwig
        File? c_a_norm_wig = if (paired_sample_m && paired_control_m) then pc.c_a_norm_wig else sc.c_a_norm_wig
        File? c_a_tdffile = if (paired_sample_m && paired_control_m) then pc.c_a_tdffile else sc.c_a_tdffile
        File? s_bigwig = if (paired_sample_m && paired_control_m) then pc.s_bigwig else sc.s_bigwig
        File? s_norm_wig = if (paired_sample_m && paired_control_m) then pc.s_norm_wig else sc.s_norm_wig
        File? s_tdffile = if (paired_sample_m && paired_control_m) then pc.s_tdffile else sc.s_tdffile
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
        Array[File?]? s_qc_statsfile = if (paired_sample_m && paired_control_m) then pc.s_qc_statsfile else sc.s_qc_statsfile
        Array[File?]? s_qc_htmlfile = if (paired_sample_m && paired_control_m) then pc.s_qc_htmlfile else sc.s_qc_htmlfile
        Array[File?]? s_qc_textfile = if (paired_sample_m && paired_control_m) then pc.s_qc_textfile else sc.s_qc_textfile
        File? s_qc_mergehtml = if (paired_sample_m && paired_control_m) then pc.s_qc_mergehtml else sc.s_qc_mergehtml
        Array[File?]? c_qc_statsfile = if (paired_sample_m && paired_control_m) then pc.c_qc_statsfile else sc.c_qc_statsfile
        Array[File?]? c_qc_htmlfile = if (paired_sample_m && paired_control_m) then pc.c_qc_htmlfile else sc.c_qc_htmlfile
        Array[File?]? c_qc_textfile = if (paired_sample_m && paired_control_m) then pc.c_qc_textfile else sc.c_qc_textfile
        File? c_qc_mergehtml = if (paired_sample_m && paired_control_m) then pc.c_qc_mergehtml else sc.c_qc_mergehtml
        Array[File?]? sp_qc_statsfile = pc.sp_qc_statsfile
        Array[File?]? sp_qc_htmlfile = pc.sp_qc_htmlfile
        Array[File?]? sp_qc_textfile = pc.sp_qc_textfile
        Array[File?]? cp_qc_statsfile = pc.cp_qc_statsfile
        Array[File?]? cp_qc_htmlfile = pc.cp_qc_htmlfile
        Array[File?]? cp_qc_textfile = pc.cp_qc_textfile
        File? s_uno_statsfile = if (paired_sample_m && paired_control_m) then pc.s_uno_statsfile else sc.s_uno_statsfile
        File? s_uno_htmlfile = if (paired_sample_m && paired_control_m) then pc.s_uno_htmlfile else sc.s_uno_htmlfile
        File? s_uno_textfile = if (paired_sample_m && paired_control_m) then pc.s_uno_textfile else sc.s_uno_textfile
        File? c_uno_statsfile = if (paired_sample_m && paired_control_m) then pc.c_uno_statsfile else sc.c_uno_statsfile
        File? c_uno_htmlfile = if (paired_sample_m && paired_control_m) then pc.c_uno_htmlfile else sc.c_uno_htmlfile
        File? c_uno_textfile = if (paired_sample_m && paired_control_m) then pc.c_uno_textfile else sc.c_uno_textfile
        File? statsfile = if (paired_sample_m && paired_control_m) then pc.statsfile else sc.statsfile
        File? htmlfile = if (paired_sample_m && paired_control_m) then pc.htmlfile else sc.htmlfile
        File? textfile = if (paired_sample_m && paired_control_m) then pc.textfile else sc.textfile
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
        File? summaryhtml = if (paired_sample_m && paired_control_m) then pc.summaryhtml else sc.summaryhtml
        File? summarytxt = if (paired_sample_m && paired_control_m) then pc.summarytxt else sc.summarytxt
        File? sp_summaryhtml = pc.sp_summaryhtml
        File? sp_summarytxt = pc.sp_summarytxt
    }
}
