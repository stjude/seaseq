#!/usr/bin/env python3
'''
========================  PEAKS ANNOTATION  ======================

USAGE
All inclusive script for Annotation of ChIP-Seq peaks.
This is included in the SEASEQ pipeline.
    ```
    python PEAKSANNO.py -g <GTF/GFF/GFF3> -c <CHROM SIZES> -p <MACS PEAKS bed> -s <MACS SUMMITS bed>
	```

The script input requirements are :
- GTF : GTF/GFF/GFF3 having either gene or transcript (prefered) annotation
- CHROM SIZES : chrom.sizes file
- MACS PEAKS bed : MACS peaks bed file
- MACS SUMMITS bed : MACS peaks summit file

Output files provided are :
	1) TSS nearest the center of peaks in __*centerofpeaks_closest.regions.txt*__
	2) peaks overlapping genes regions in __*peaks_within_genebody.regions.txt*__
	3) peaks overlapping promoters in __*peaks_within_promoter.regions.txt*__
	4) peaks overlapping windows in __*peaks_within_window.regions.txt*__
	5) peaks identified in previous overlapping regions and comparison of
           all regions in __*peaks_compared_regions.peaks.txt*__
	6) genes identified in previous overlapping regions and comparison of
           all regions in __*peaks_compared_regions.genes.txt*__
	7) distribution graphs in __*peaks_compared_regions.distribution.pdf*__

Definitions of annotation regions:
	- TSS: TSS (transcription start site)
	- promoters: 1kb +/- TSS
	- window: 10kb upstream --> 3kb downstream of genebody
	- genebody: 1kb upstream of TSS --> TES (transcription end site)

Dependencies of script:
	* bedtools
	* python3

==================================================================
'''

import os
import sys
import re
import argparse
import subprocess
import matplotlib.pyplot as plt
from pathlib import Path
import collections
from collections import defaultdict

if not os.path.exists('annotation'):
    os.makedirs('annotation')

def parse_genelocations(chromz, results, up_dist, down_dist, tss=True):
    """ Parse genomic regions
    Args:
        chromz (dict) : Chromosomal sizes
        results (string) : Individual gene coordinates
        up_dist (int) : genomic distance upstream / +
        down_dist (int) : genomic distance downstream / -
        tss (boolean) : TSS or entire gene (true if tss)
    """

    lines = results.split("\t")
    lines[3] = int(lines[3])
    lines[4] = int(lines[4])

    if lines[6] == "+":
        if tss:
            end_region = 4
        else:
            end_region = 3

        end_position = lines[end_region] + down_dist
        start_position = lines[3] - up_dist

    elif lines[6] == "-":
        if tss:
            end_region = 3
        else:
            end_region = 4

        start_position = lines[end_region] - down_dist
        end_position = lines[4] + up_dist

    if start_position < 1:
        start_position = 1

    if end_position > int(chromz[lines[0]]):
        end_position = chromz[lines[0]]

    output_gff = ("{0}\t{1}\t{2}\t{3}\n".format("\t".join(lines[0:3]),
                                                start_position, end_position, "\t".join(lines[5:])))
    return output_gff


def gtf_to_genes(chrom_sizes_file, gtf_file):
    """ GTF to GENES
    Args:
        chrom_sizes_file (file) : Chromosomal sizes file
        gtf_file (file) : GTF file
    """
    #initialize outputfiles
    original_output = open('annotation/original_genes.gff', 'w')
    promoters_output = open('annotation/promoter-region_genes.gff', 'w')
    window_output = open('annotation/window-region_genes.gff', 'w')
    tss_output = open('annotation/TSS-region_genes.gff', 'w')
    genebody_output = open('annotation/genebody-region_genes.gff', 'w')

    haschr = False #check 'chr' prefix in gtf file
    
    #read ucsc chrom sizes, gtf
    chrom_sizes_input = open(chrom_sizes_file, 'r')
    chrom_sizes_dict = {}
    for line in chrom_sizes_input:
        lines = line.split('\t')
        chrom_sizes_dict[lines[0]] = lines[1].rstrip("\n")
        if lines[0].startswith('chr'):
            haschr = True

    #feature = options.feature
    feature_dict = {}

    #reading the feature type to be used
    gtf_input = open(gtf_file, 'r')
    for line in gtf_input:
        if not line.startswith('#'):
            lines = line.split("\t")
            feature_dict[lines[2]] = lines[2]

    if 'transcript' in feature_dict:
        feature = "transcript"
    elif 'gene' in feature_dict:
        feature = "gene"
    else:
        sys.exit("ERROR :\tGTF/GFF with either transcript/gene annotation is needed")
    print("NOTE :\tFeature type used is '%s'" %(feature))

    #organize gtf file
    gtf_input = open(gtf_file, 'r')
    for line in gtf_input:
        if not line.startswith('#'):
            lines = line.rstrip("\n").split("\t")
            if haschr and not lines[0].startswith('chr'):
                lines[0] = "chr"+lines[0]
            elif not haschr and lines[0].startswith('chr'):
                lines[0] = lines[0][3:]
            if lines[2] == feature:
                newline = lines[8].split(';')
                
                if gtf_file.split('.')[-1] == 'gff' or gtf_file.split('.')[-1] == 'gff3':
                    if re.search(";transcript_id", lines[8]):
                        transcript_id = [s for s in newline if "transcript_id=" in s][0]
                    else:
                        transcript_id = newline[0]
                    if re.search(";gene_name", lines[8]):
                        gene_name = [s for s in newline if "gene_name=" in s][0]
                    else:
                        gene_name = newline[0]
                    transcript_id = re.sub('[\"\;]', '', transcript_id.split('=')[-1]) #clean transcript_id
                    gene_name = re.sub('[\"\;]', '', gene_name.split('=')[-1]) #clean gene_name
                elif gtf_file.split('.')[-1] == 'gtf':
                    if re.search("; transcript_id", lines[8]):
                        transcript_id = [s for s in newline if "transcript_id " in s][0]
                    else:
                        transcript_id = newline[0]
                    if re.search("; gene_name", lines[8]):
                        gene_name = [s for s in newline if "gene_name " in s][0]
                    else:
                        gene_name = newline[0]
                    transcript_id = re.sub('[\"\;]', '', transcript_id.split(' ')[-1]) #clean transcript_id
                    gene_name = re.sub('[\"\;]', '', gene_name.split(' ')[-1]) #clean gene_name
                    
                results = ("{0}\t{1}\t{2}\t{3}".format(lines[0],
                                                       "\t".join(lines[1:8]),
                                                       gene_name, transcript_id))

                try:
                    if chrom_sizes_dict[lines[0]]:
                        original_output.write(results + "\n")
                        final_output = parse_genelocations(chrom_sizes_dict, results, 1000, 1000, False)
                        promoters_output.write(final_output)
                        final_output = parse_genelocations(chrom_sizes_dict, results, 10000, 3000, True)
                        window_output.write(final_output)
                        final_output = parse_genelocations(chrom_sizes_dict, results, 0, 0, False)
                        tss_output.write(final_output)
                        final_output = parse_genelocations(chrom_sizes_dict, results, 1000, 0, True)
                        genebody_output.write(final_output)
                except KeyError:
                    continue


def annotate_regions(region_file, output_file, macs_summits_file=None):
    """ Annotate Peaks based on genomic features
    Args:
        region_file (file) : Region file
        output_file (string) : Output file name
        macs_summits_file (file) : MACS summits file
    """

    #initialize chromosomal positions based on gff provided
    opt_chrom = 0
    opt_start = 1
    opt_stop = 2
    opt_peakid = 3
    opt_score = 4
    opt_none = 6
    opt_gene = 8
    opt_transcript = 9

    summits_dict = {}
    if macs_summits_file:
        #read summits file
        summits_input = open(macs_summits_file, 'r')
        for line in summits_input:
            line = line.strip('\n').split('\t')
            string_identifier = line[opt_start], line[opt_stop]
            summits_dict[line[opt_peakid]] = string_identifier

    #read intersectbed file
    intersect_input = open(region_file, 'r')
    peaks_dict = {}
    peaks_score = {}
    gene_dict = {}
    transcript_dict = {}
    peak_region = []
    for line in intersect_input:
        line = line.strip('\n').split('\t')
        #print(line[opt_none],"\n")
        if int(line[opt_none]) != -1:
            string_identifier = line[opt_chrom], line[opt_start], line[opt_stop]
            peaks_dict[string_identifier] = line[opt_peakid]
            peaks_score[line[opt_peakid]] = line[opt_score]
            peak_region.append(string_identifier)

            if string_identifier in gene_dict:
                gene_dict[string_identifier] = ("{0},{1}".format(gene_dict[string_identifier],
                                                                 line[opt_gene]))
                transcript_dict[string_identifier] = ("{0},{1}".
                                                      format(transcript_dict[string_identifier],
                                                             line[opt_transcript]))
            else:
                gene_dict[string_identifier] = line[opt_gene]
                transcript_dict[string_identifier] = line[opt_transcript]

    #each peak id is unique
    unique_peak_region = []
    for  eachpeak_region in peak_region:
        if eachpeak_region not in unique_peak_region:
            unique_peak_region.append(eachpeak_region)

    #write to output file
    final_output = open(output_file, 'w')

    if macs_summits_file:
        final_output.write("peak_id\tchrom\tpeak_start\tpeak_end\t"
                           "summit_start\tsummit_end\tgene_name\t"
                           "transcript_id\tpeak_score\n")

        #print regions identified
        for chromregion in unique_peak_region:
            peakid = peaks_dict[chromregion]
            final_output.write("%s\t%s\t%s\t%s\t%s\t%s\n"
                               %(peakid, '\t'.join(chromregion), '\t'.join(summits_dict[peakid]),
                                 ','.join(unique(gene_dict[chromregion])),
                                 ','.join(unique(transcript_dict[chromregion])),
                                 peaks_score[peakid]))
    else:
        final_output.write("peak_id\tchrom\tpeak_start\tpeak_end\t"
                           "gene_name\ttranscript_id\tpeak_score\n")

        #print regions identified
        for chromregion in unique_peak_region:
            peakid = peaks_dict[chromregion]
            final_output.write("%s\t%s\t%s\t%s\t%s\n"
                               %(peakid, '\t'.join(chromregion),
                                 ','.join(unique(gene_dict[chromregion])),
                                 ','.join(unique(transcript_dict[chromregion])),
                                 peaks_score[peakid]))

def unique(list1):
    """ Get unique list
    Args:
        list1 (list) : list of variable separated by comma ','

    """

    # intilize a null list
    unique_list = []

    # traverse for all elements
    for variable in list1.split(","):
        if variable not in unique_list:
            unique_list.append(variable)
    return unique_list


def include_genes(macs_peaks_file, file_name):
    """ Annotate Genes in peak regions
    Args:
        macs_peaks_file (file) : MACS peaks file
        output_file (string) : Output file name
    """
    peaks_input = open(macs_peaks_file, 'r')
    peaks_chr_dict = {}
    peaks_count = 0
    for line in peaks_input:
        line = line.strip('\n').split('\t')
        peaks_count += 1
        peaks_chr_dict[line[3]] = (line[0:3])

    #reading the gtf genes file
    gene_names_dict = {}
    gene_names_input = open("tempPrefix.genes.names", 'r')
    for line in gene_names_input:
        line = line.strip('\n').split('\t')
        gene_names_dict[(line[0], line[1])] = line[0]

    #reading the genomic regions files
    region_names_list = []
    region_dict = defaultdict(dict)
    peaks_id_dict = defaultdict(dict)
    for region_prefix in ("promoter", "genebody", "window", "closest"):
        region = region_prefix + ".tempPrefix"
        region_name = Path(region).stem
        region_names_list.append(region_name.upper())
        region_input = open(region, 'r')
        for line in region_input:
            line = line.strip('\n').split('\t')
            region_dict[region_name.upper()][(line[8], line[9])] = (line[8], line[9])
            peaks_id_dict[region_name.upper()][line[3]] = (line[0:3])
    peaks_id_list = peaks_chr_dict.keys()
    peaks_id_list = sorted(peaks_id_list, key=lambda item:
                           (int(item[10:]) if item[10:].isdigit()
                            else print(item[10:], type(item[10:])), item))

    #order regions and gene names
    region_names_list.sort()
    orderedgene_names_dict = collections.OrderedDict(sorted(gene_names_dict.items()))

    #barcoding the identified regions
    barcoderegions_dict = {}
    master_generegion_dict = defaultdict(dict)
    for eachgene in orderedgene_names_dict:
        for eachkey in region_names_list:
            found = 0
            if eachgene in region_dict[eachkey]:
                found = 1
            master_generegion_dict[eachkey][eachgene] = found
            if eachgene in barcoderegions_dict:
                barcoderegions_dict[eachgene] += str(found)
            else:
                barcoderegions_dict[eachgene] = str(found)

    barcodepeaks_id_dict = {}
    master_peaksid_dict = defaultdict(dict)
    for eachpeaks_id in peaks_id_list:
        for eachkey in region_names_list:
            found = 0
            #print(eachpeaks_id)
            if eachpeaks_id in peaks_id_dict[eachkey]:
                found = 1
                #print("yes")
            master_peaksid_dict[eachkey][eachpeaks_id] = found
            if eachpeaks_id in barcodepeaks_id_dict:
                barcodepeaks_id_dict[eachpeaks_id] += str(found)
            else:
                barcodepeaks_id_dict[eachpeaks_id] = str(found)

    #write to output file
    gene_file_name = file_name.replace('.txt', '.genes.txt')

    final_output = open(gene_file_name, 'w')
    final_output.write("GeneName\tTranscript_ID\tBARCODE\t%s\tEVER\n"
                       %('\t'.join(region_names_list)))

    for eachgene in orderedgene_names_dict:
        string_identifier = [barcoderegions_dict[eachgene]]
        found_sum = 0
        for eachkey in region_names_list:
            string_identifier.append(str(master_generegion_dict[eachkey][eachgene]))
            found_sum += master_generegion_dict[eachkey][eachgene]
        if found_sum >= 1:
            string_identifier.append("1")
        else:
            string_identifier.append("0")

        final_output.write("%s\t%s\n" %('\t'.join(eachgene), '\t'.join(string_identifier)))

    peaks_file_name = file_name.replace('.txt', '.peaks.txt')
    pdf_file_name = file_name.replace('.txt', '.distribution.pdf')
    final_output_peaks = open(peaks_file_name, 'w')
    final_output_peaks.write("chrom\tstart\tstop\tpeak_id\tBARCODE\t%s\n"
                             %('\t'.join(region_names_list)))
    for eachpeaks_id in peaks_id_list:
        string_identifier = [barcodepeaks_id_dict[eachpeaks_id]]
        for eachkey in region_names_list:
            string_identifier.append(str(master_peaksid_dict[eachkey][eachpeaks_id]))

        final_output_peaks.write("%s\t%s\t%s\n"
                                 %('\t'.join(peaks_chr_dict[eachpeaks_id]),
                                   eachpeaks_id, '\t'.join(string_identifier)))

    #the bar plot
    region_count_dict = {}
    #for eachregion in region_count_dict:
    region_count_dict["Promoter\n(1kb up/down TSS)"] = round(
        len(peaks_id_dict["PROMOTER"])*100/peaks_count)
    region_count_dict["Genebody\n(1kb up TSS/TES)"] = round(
        len(peaks_id_dict["GENEBODY"])*100/peaks_count)
    region_count_dict["Window\n(10kb up TSS/3kb down TES)"] = round(
        len(peaks_id_dict["WINDOW"])*100/peaks_count)

    rects = plt.bar(region_count_dict.keys(), region_count_dict.values(), width=0.4)
    axes = plt.gca()
    axes.set_ylim([0, 100])
    peaks_title = Path(macs_peaks_file).stem
    plt.title(peaks_title.split('.sorted')[0])
    plt.ylabel("% of peaks")

    #label the bars
    xpos = 'center'
    ha = {'center': 'center', 'right': 'left', 'left': 'right'}
    offset = {'center': 0.5, 'right': 0.57, 'left': 0.43}  # x_txt = x + w*off

    for rect in rects:
        height = rect.get_height()
        plt.text(rect.get_x() + rect.get_width()*offset[xpos], 1.01*height,
                 '{}%'.format(height), ha=ha[xpos], va='bottom')

    #save figure
    plt.savefig(pdf_file_name)

    for eachregion in ("PROMOTER", "GENEBODY", "WINDOW"):
        print("Total number of %s = %d" %(eachregion, len(peaks_id_dict[eachregion])))
        print("Percentage of %s = %.3f" %(eachregion,
                                          len(peaks_id_dict[eachregion])*100/peaks_count))


def main():
    '''
    Generate genomic coordinates of all promoters, window, genebody, TSS

    Definitions: promoters: +/- 1kb of TSS (transcription start site)
                 window: 10kb upstream --> 3kb downstream of genebody
                 genebody: 1kb upstream of TSS --> TTS/TES (transcription termination/end site)
                 TSS: TSS (transcription start site)

    Annotate Peaks based on genomic features listed above

    Comparison master files for all genes/transcripts
        and their genomic regions within peaks called.
    Genomic regions are specified based on regions provided,
        default expectation are promoter, genebody, window, closest

    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gtf", dest="gtf", required=True,
                        help="Enter .gtf/gff file to be processed.")
    parser.add_argument("-c", "--chrom", dest="chrom_sizes", required=True,
                        help="Enter ucsc chrom sizes file to be processed.")
    parser.add_argument("-p", "--peaks", dest="macs_peaks", required=True,
                        help="Enter MACS peaks bed file.")
    parser.add_argument("-s", "--summit", dest="macs_summits", required=False,
                        help="Enter MACS summit bed file.")

    options = parser.parse_args()
    #print(options)

    # Stage 1
    # Extract genomic regions of interest
    gtf_to_genes(options.chrom_sizes, options.gtf)

    # Turn 4 column bed to 5 (for SICER peaks: include unique id)
    command = "wc -l " + options.macs_peaks + " | awk -F' '  '{print $1}'"
    checkcolumns = int(subprocess.check_output(command,shell=True).strip())
    if checkcolumns > 0:
        command = "awk -F'\\t' '{print NF; exit}' " + options.macs_peaks
        numberofcolumns = int(subprocess.check_output(command,shell=True).strip())
    else:
        sys.exit("ERROR :\tNo peaks were identified\n")

    command = "head -n 1 " + options.macs_peaks + ' | cut -f5 '
    fifthcol_value = (str(subprocess.check_output(command,shell=True).strip()).split("'"))[1] #FDRisland from SICER produce a fifth empty column
    
    new_macs_peaks = options.macs_peaks

    if numberofcolumns == 4 or len(fifthcol_value) < 1:
        new_macs_peaks = options.macs_peaks + ".tempPrefix"
        # command = "cat " + options.macs_peaks + ' | awk -F\\\\t ' + "'" + \
        #           '{print $1 "\\t" $2 "\\t" $3 "\\t" $1":"$2"-"$3 "\\t" $4}' + \
        #           "' > tempPrefix-new.macs_peaks.bed"
        # print(command)
        # os.system(command)
        macs_peaks_output = open(new_macs_peaks, 'w')
        macs_peaks_input = open(options.macs_peaks,'r')
        macs_peaks_count = 1
        for line in macs_peaks_input:
            line = line.strip('\n').split('\t')
            new_peaks_id = "SICERpeak_" + str(macs_peaks_count)
            results = ("{0}\t{1}\t{2}\n".format("\t".join(line[0:3]), new_peaks_id, line[3]))
            macs_peaks_output.write(results)
            macs_peaks_count += 1
        macs_peaks_output.close()
        
    # Convert TSS to compatible bedtools (closestBed) format and obtain center of peaks coordinate
    command = "cut -f1,4,5,9,10 annotation/TSS-region_genes.gff | sort -k1,1 -k2,2n " + \
              "> tempPrefix-sorted.TSS; " + \
              "cat " + new_macs_peaks + ' | awk -F\\\\t ' + "'" + \
              '{print $1 "\\t" int(($2+$3)/2) "\\t" int(($2+$3)/2) "\\t" $4 "\\t" $5}' + \
              "' | sort -k1,1 -k2,2n > tempPrefix-centerofpeaks.bed;" + \
              "closestBed -t all -a tempPrefix-centerofpeaks.bed -b tempPrefix-sorted.TSS" + \
              "> closest.tempPrefix"
    os.system(command)

    # Stage 2
    # Annotate peaks based on previously extracted genomic coordinates
    for region in ("promoter", "genebody", "window"):
        command = 'intersectBed -b ' + new_macs_peaks + ' -a annotation/' + region + \
                  '-region_genes.gff -wb | cut -f1,4,5,9,10,11,12,13,14,15 |' + "awk -F\\\\t '" + \
                  '{print $6 "\\t" $7 "\\t" $8 "\\t" $9 "\\t" $10 "\\t" $1 "\\t" $2 "\\t" $3 ' + \
                  '"\\t" $4 "\\t" $5}' + \
                  "' > " + region + ".tempPrefix"
        os.system(command)

    # Convert annotated regions to tab delimited file and include summit location
    for region in ("promoter", "genebody", "window", "closest"):
        region_file = region + ".tempPrefix"
        output_file = "peaks_within_" + region + ".regions.txt"
        if options.macs_summits:
            annotate_regions(region_file, output_file, options.macs_summits)
        else:
            annotate_regions(region_file, output_file)

    command = "mv peaks_within_closest.regions.txt centerofpeaks_closest.regions.txt; \
              cut -f 9,10 annotation/original_genes.gff > tempPrefix.genes.names"
    os.system(command)

    # Stage 3
    # Binary matrix of all genes
    include_genes(new_macs_peaks, "peaks_compared_regions.txt")

    command = "rm -rf *tempPrefix*"
    #os.system(command)


if __name__ == "__main__":
    main()
