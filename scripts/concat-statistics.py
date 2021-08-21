#!/usr/bin/env python3
#Concatenate stats of all samples provided

import argparse

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sample", dest="sample", required=True,
                        help="Enter sample seaseq config ml.")
    parser.add_argument("-c", "--control", dest="control", required=True,
                        help="Enter control seaseq config ml.")
    parser.add_argument("-q", "--overall", dest="overall", required=True,
                        help="Enter overall final seaseqconfig ml.")
    parser.add_argument("-o", "--outfile", dest="outfile", required=False,
                        default="summarystats-stats.htmlx",
                        help="Enter output file name.")

    options = parser.parse_args()

    #all statistics
    stats = {
        1 : 'Overall_Quality',
        2 : 'Raw_Reads',
        3 : 'Base_Quality',
        4 : 'Sequence_Diversity',
        5 : 'Aligned_Percent',
        6 : 'Estimated_Fragment_Width',
        7 : 'Estimated_Tag_Length',
        8 : 'NRF',
        9 : 'PBC',
        10 : 'NSC',
        11 : 'RSC',
        12 : 'FRiP',
        13 : 'Total_Peaks',
        14 : 'Normalized_Peaks*',
        15 : 'Linear_Stitched_Peaks',
        16 : 'SE-like_Enriched_Regions'
    }

    #color scheme
    color = {
        -2 : "#FF0000",
        -1 : "#FF8C00",
        0 : "#FFFF00",
        1 : "#ADFF2F",
        2 : "#008000"
    }

    #merged statistics
    mergestats = [1, 14, 15, 16]
    
    #initialize dictionaries
    SQCvalue = {}
    CQCvalue = {}
    OQCvalue = {}
    SQCscore = {}
    CQCscore = {}
    OQCscore = {}

    if options.outfile.endswith("stats.htmlx"):
        statsfile = options.outfile
    else:
        statsfile = options.outfile + '-stats.htmlx'

    htmlfile = statsfile
    textfile = htmlfile.replace("stats.htmlx", "stats.txt")

    sample_content = open(options.sample, "r")
    control_content = open(options.control, "r")
    overall_content = open(options.overall, "r")

    for line in sample_content:
        eachdata = line.rstrip('\n').split('\t')
        SQCvalue[eachdata[0]] = eachdata[1]
        SQCscore[eachdata[0]] = int(eachdata[2])
    for line in control_content:
        eachdata = line.rstrip('\n').split('\t')
        CQCvalue[eachdata[0]] = eachdata[1]
        CQCscore[eachdata[0]] = int(eachdata[2])
    for line in overall_content:
        eachdata = line.rstrip('\n').split('\t')
        OQCvalue[eachdata[0]] = eachdata[1]
        OQCscore[eachdata[0]] = int(eachdata[2])

    htmlheader = "<table class='results'><tr><th>DATA"
    textheader = "DATA"
    samplehtmlvalues = "<tr><td><center>SAMPLE</center></td>"
    controlhtmlvalues = "<tr><td><center>CONTROL</center></td>"
    sampletextvalues = "SAMPLE"
    controltextvalues = "CONTROL"

    for key in sorted(stats.keys()):
        #change space to underscore for txt fileprint(key)
        convertheader = stats[key].replace('_', ' ')
        textheader += "\t" + stats[key]
        htmlheader += "</th><th>" + convertheader
        if key in mergestats:
            if key == 14:
                samplehtmlvalues += "<td rowspan='2' bgcolor='" + \
                                    color[OQCscore['Total_Peaks']] + "'><center>" + \
                                    OQCvalue['Total_Peaks'] + "</center></td>"
                sampletextvalues += "\t" + OQCvalue['Total_Peaks']
                controltextvalues += "\t"
            else:
                samplehtmlvalues += "<td rowspan='2' bgcolor='" + \
                                    color[OQCscore[stats[key]]] + "'><center>" + \
                                    OQCvalue[stats[key]] + "</center></td>"
                sampletextvalues += "\t" + OQCvalue[stats[key]]
                controltextvalues += "\t"
        else:
            samplehtmlvalues += "<td bgcolor='" + color[SQCscore[stats[key]]] + \
                                "'><center>" + SQCvalue[stats[key]] + "</center></td>"
            controlhtmlvalues += "<td bgcolor='" + color[CQCscore[stats[key]]] + \
                                "'><center>" + CQCvalue[stats[key]] + "</center></td>"
            sampletextvalues += "\t" + SQCvalue[stats[key]]
            controltextvalues += "\t" + CQCvalue[stats[key]]

    htmlheader += "</th></tr>"
    samplehtmlvalues += "</tr>"
    controlhtmlvalues += "</tr>"

    writehtmlfile = open(htmlfile, 'w')
    writetextfile = open(textfile, 'w')

    writehtmlfile.write(htmlheader + "\n" + samplehtmlvalues + "\n" + controlhtmlvalues)
    writetextfile.write(textheader + "\n" + sampletextvalues + "\n" + controltextvalues)

if __name__ == "__main__":
    main()
