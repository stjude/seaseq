version 1.0

import "../tasks/bedtools.wdl"
import "../tasks/bowtie.wdl"
import "../tasks/samtools.wdl"

workflow mapping {
    input {
        File fastqfile
        File? fastqfile_R2
        Int insert_size = 600
        Array[File] index_files
        File? metricsfile
        File? blacklist
        Int? read_length = 75
        Boolean paired_end = false
        String? strandedness = 'fr'
        String default_location = "BAM_files"
        String? results_name
    }

    if ( defined(metricsfile) ) {
        call bowtie.bowtie as SE_map {
            input :
                fastqfile=fastqfile,
                index_files=index_files,
                metricsfile=metricsfile,
                read_length=read_length
        }
    }

    if ( defined(fastqfile_R2) ) {
        call bowtie.bowtie as PE_map {
            input :
                fastqfile=fastqfile,
                fastqfile_R2=fastqfile_R2,
                index_files=index_files,
                insert_size=insert_size,
                strandedness=strandedness,
                prefix=results_name
        }
    }   

    call samtools.viewsort {
        input :
            samfile=select_first([SE_map.samfile, PE_map.samfile]),
            default_location=default_location,
            paired_end=paired_end
    }

    call samtools.indexstats {
        input :
            bamfile=viewsort.sortedbam,
            default_location=default_location
    }

    if ( defined(blacklist) ) {
        # remove blacklist regions
        String string_blacklist = "" #buffer to allow for blacklist optionality
        File blacklist_file = select_first([blacklist, string_blacklist])
        if ( defined(metricsfile) ) {
            call bedtools.intersect as rmblklist {
                input :
                    fileA=viewsort.sortedbam,
                    fileB=blacklist_file,
                    default_location=default_location,
                    nooverlap=true
            }
        }

        if ( defined(fastqfile_R2) ) {
            String string_fixmate = "" #buffer to allow for blacklist optionality
            File fixmate_file = select_first([viewsort.fixmatebam, string_fixmate])
            call bedtools.pairtobed as pairtobed {
                input :
                    fileA=fixmate_file,
                    fileB=blacklist_file,
                    default_location=default_location
            }
        }

        File bklist_bamfile = select_first([rmblklist.intersect_out, pairtobed.pairtobed_out])

        call samtools.indexstats as bklist {
            input :
                bamfile=bklist_bamfile,
                default_location=default_location
        }          
    } # end if (blacklist provided)

    call samtools.markdup {
        input :
            bamfile=select_first([bklist_bamfile, viewsort.sortedbam]),
            default_location=default_location
    }

    call samtools.indexstats as mkdup {
        input :
            bamfile=markdup.mkdupbam,
            default_location=default_location
    }

    File parse_sorted_bam = select_first([viewsort.fixmatebam, viewsort.sortedbam])
    
    output {
        File sorted_bam = viewsort.sortedbam
        File as_sortedbam = parse_sorted_bam
        File bam_index = indexstats.indexbam
        File bam_stats = indexstats.flagstats
        File? bklist_bam = bklist_bamfile
        File? bklist_index = bklist.indexbam
        File? bklist_stats = bklist.flagstats
        File mkdup_bam = markdup.mkdupbam
        File mkdup_index = mkdup.indexbam
        File mkdup_stats = mkdup.flagstats
    }
}
