version 1.0

import "https://raw.githubusercontent.com/adthrasher/seaseq/refactor/workflows/tasks/bedtools.wdl"
import "https://raw.githubusercontent.com/adthrasher/seaseq/refactor/workflows/tasks/bowtie.wdl"
import "https://raw.githubusercontent.com/adthrasher/seaseq/refactor/workflows/tasks/samtools.wdl"

workflow mapping {
    input {
        File fastqfile
        Array[File] index_files
        File metricsfile
        File? blacklist
        Int? read_length = 75
        String default_location = "BAM_files"
    }

    call bowtie.bowtie {
        input :
            fastqfile=fastqfile,
            index_files=index_files,
            metricsfile=metricsfile,
            read_length=read_length
    }   

    call samtools.viewsort {
        input :
            samfile=bowtie.samfile,
            default_location=default_location
    }

    call samtools.indexstats {
        input :
            bamfile=viewsort.sortedbam,
            default_location=default_location
    }

    if ( defined(blacklist) ) {
        # remove blacklist regions
        String string_blacklist = "" #buffer to allow for blacklist optionality
        File blacklist_ = select_first([blacklist, string_blacklist])
        call bedtools.intersect as rmblklist {
            input :
                fileA=viewsort.sortedbam,
                fileB=blacklist_,
                default_location=default_location,
                nooverlap=true
        }
        call samtools.indexstats as bklist {
            input :
                bamfile=rmblklist.intersect_out,
                default_location=default_location
        }
    } # end if (blacklist provided)

    call samtools.markdup {
        input :
            bamfile=select_first([rmblklist.intersect_out, viewsort.sortedbam]),
            default_location=default_location
    }

    call samtools.indexstats as mkdup {
        input :
            bamfile=markdup.mkdupbam,
            default_location=default_location
    }

    output {
        File sorted_bam = viewsort.sortedbam
        File bam_index = indexstats.indexbam
        File bam_stats = indexstats.flagstats
        File? bklist_bam = rmblklist.intersect_out
        File? bklist_index = bklist.indexbam
        File? bklist_stats = bklist.flagstats
        File mkdup_bam = markdup.mkdupbam
        File mkdup_index = mkdup.indexbam
        File mkdup_stats = mkdup.flagstats
    }
}
