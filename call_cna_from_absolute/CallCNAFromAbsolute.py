# =====================================
# CallCNAFromAbsolute.py
# Version 1.0
# David Liu, Van Allen Lab
# davidliu@broadinstitute.org
# =====================================

# Given gene annotated segtab files from ABSOLUTE generate amplifications, LOH, and deletions
# based on process from Brastianos, Carter et al Cancer Discovery 2015

import pandas as pd

import os
import numpy as np
import glob
import csv
import sys
import math
import argparse
from collections import defaultdict

def GenCNADistribution(geneSegDF):
    """iterates through pandas dataframe of the geneSegFile and generates a dictionary of copy number keys and
    frequency values"""
    histCNAs = {}
    segConsidered = {} # to keep track of segments we have already considered
    totalBP = 0
    totalSeg = 0

    for i in range(0, len(geneSegDF)):
        chromosome = str(geneSegDF.ix[i, 'chr'])
        segStart = geneSegDF.ix[i, 'start']
        segEnd = geneSegDF.ix[i, 'end']
        rCNA1 = geneSegDF.ix[i, 'rescaled.cn.a1']
        rCNA2 = geneSegDF.ix[i, 'rescaled.cn.a2']

        if chromosome in segConsidered:
            if (str(segStart) + "-" + str(segEnd)) in segConsidered[chromosome]: # already previously considered
                continue # skip to next entry
            else: # add to list of processed segments
                segConsidered[chromosome][str(segStart) + "-" + str(segEnd)] = 1
        else: # first time seeing the chromosome!  
            segConsidered[chromosome] = {}
            segConsidered[chromosome][str(segStart) + "-" + str(segEnd)] = 1

        # Consider special case where rescaled.CN is "NA" -- this is case where the number of probes is too small. 
        # To preserve events, will assume uneven distributions; rCNA1 = 1 if total CN >=2 and 0 otherwise.
        if pd.isnull(rCNA1):
            expectedTotal = geneSegDF.ix[i, 'expected_total_cn']
            if expectedTotal >= 2:
                rCNA1 = 1.0
            else:
                rCNA1 = 0.0
            rCNA2 = expectedTotal - rCNA1

        if str(rCNA1) in histCNAs: # add # of base pairs
            histCNAs[str(rCNA1)]['bp'] += segEnd - segStart + 1
        else: # not seen previously
            histCNAs[str(rCNA1)] = {}
            histCNAs[str(rCNA1)]['bp'] = segEnd - segStart + 1

        if str(rCNA2) in histCNAs: # add # of base pairs
            histCNAs[str(rCNA2)]['bp'] += segEnd - segStart + 1
        else: # not seen previously
            histCNAs[str(rCNA2)] = {}
            histCNAs[str(rCNA2)]['bp'] = segEnd - segStart + 1

        totalBP += 2*(segEnd - segStart + 1)
        totalSeg += 1

    # Now iterate through the histogram to generate percentiles at each copy number
    cns = sorted([float(i) for i in histCNAs.keys()])

    fractionBelow = 0
    for i in range(0, len(cns)):
        cn = cns[i]
        histCNAs[str(cn)]['fractionBelow'] = fractionBelow
        histCNAs[str(cn)]['fraction'] = histCNAs[str(cn)]['bp']/float(totalBP)
        sys.stdout.write("loading copy number {}: {} {}\n".format(cn, fractionBelow, histCNAs[str(cn)]['bp']/float(totalBP)))
        fractionBelow += histCNAs[str(cn)]['bp']/float(totalBP)

    sys.stdout.write("fractionBelow should equal 1: {}\n".format(fractionBelow))
    sys.stdout.write("Total number of segments considered: {}\n".format(totalSeg))
    sys.stdout.write("Total number of basepairs: {}\n".format(totalBP))
    sys.stdout.write(str(segConsidered['1']))
    
    return histCNAs


def CalcFocality(cn, CNAHist):
    """Given a copy number and a distribution of copy numbers, determines the fraction of the genome that is
     < or > the value (whichever is smaller), returns (fractionbelow, fractionabove, focality)"""
    if str(cn) not in CNAHist:
        sys.stdout.write("{} not found in copy number histogram.".format(cn))
        exit()

    fractionBelow = CNAHist[str(cn)]['fractionBelow'] 
    fraction = CNAHist[str(cn)]['fraction']
    fractionAbove = 1 - fractionBelow - fraction

    if fractionBelow < fractionAbove:
        return (fractionBelow, fractionAbove, 1-fractionBelow)
    else:
        return (fractionBelow, fractionAbove, 1-fractionAbove)


def GenFocality(geneSegDF, CNAHist):
    """Given a df of genes and rCN, calculate the focality = fraction of genome that has < (deletion) or > (amp)
    the given copy number. Then makes the CNA call, and Saves it back as additional columns in DF."""
    fractionBelow1 = []
    fractionBelow2 = []
    fractionAbove1 = []
    fractionAbove2 = []

    focality1 = []
    focality2 = []

    call1 = []
    call2 = []

    for i in range(0, len(geneSegDF)):
        rCNA1 = geneSegDF.ix[i, 'rescaled.cn.a1']
        rCNA2 = geneSegDF.ix[i, 'rescaled.cn.a2']
        
        # Special case where rescaled.CN is nan due to low probe number.
        if pd.isnull(rCNA1): # We want to maintain focal events; if expected CN > 2 then we will make rCN1 1 otherwise 0 and rCN2 the remainder
            expectedTotal = geneSegDF.ix[i, 'expected_total_cn']
            if expectedTotal >= 2:
                rCNA1 = 1.0
            else: 
                rCNA1 = 0.0
            rCNA2 = expectedTotal - rCNA1
        
        (below1, above1, focal1) = CalcFocality(rCNA1, CNAHist)
        (below2, above2, focal2) = CalcFocality(rCNA2, CNAHist)

        focality1.append(focal1)
        focality2.append(focal2)

        fractionBelow1.append(below1)
        fractionBelow2.append(below2)

        fractionAbove1.append(above1)
        fractionAbove2.append(above2)

        call1.append(GenCNACall(rCNA1, focal1))
        call2.append(GenCNACall(rCNA2, focal2))

    # now add the columns to the dataframe

    geneSegDF['focality_1'] = pd.Series(focality1, index=geneSegDF.index)
    geneSegDF['focality_2'] = pd.Series(focality2, index=geneSegDF.index)

    geneSegDF['fr_below_1'] = pd.Series(fractionBelow1, index=geneSegDF.index)
    geneSegDF['fr_below_2'] = pd.Series(fractionBelow2, index=geneSegDF.index)

    geneSegDF['fr_above_1'] = pd.Series(fractionAbove1, index=geneSegDF.index)
    geneSegDF['fr_above_2'] = pd.Series(fractionAbove2, index=geneSegDF.index)

    geneSegDF['called_CNA1'] = pd.Series(call1, index=geneSegDF.index)
    geneSegDF['called_CNA2'] = pd.Series(call2, index=geneSegDF.index)


def GenCNACall(rCN, focality):
    """given rCN and focality generate the alteration"""
    # deletion
    if rCN < 0.25 and focality > 0.995:
        return "del"

    # high level amp
    if focality > (0.98 - 1/float(7) * math.log(rCN/float(7),2)):
        return "high amp"
    
    # amplification
    if focality > (0.98 - 0.2 * math.log(rCN/float(5), 2)):
        return "amp"

    return "none"


def OutputArmSummary(geneSegDF, outputDir, filename):
    """Summarize which arms and arm segments exhibit amplification, high amplification, or deletion"""
    summaryDict = {'amp': set(), 'high amp': set(), 'del': set()}
    for index, row in geneSegDF.iterrows():
        chrom = row['chr']
        effect_1 = row['called_CNA1']
        effect_2 = row['called_CNA2']
        band = row['band']
        arm = row['arm']
        if effect_1 is not 'none':
            summaryDict[effect_1].add('{}{}\t{}{}'.format(chrom, arm, chrom, band))
        if effect_2 is not 'none':
            summaryDict[effect_2].add('{}{}\t{}{}'.format(chrom, arm, chrom, band))
    of = open(os.path.join(outputDir, '{}.cna_processed_arm_summary.tsv'.format(filename)), 'w')
    for effect, chrom_sections in summaryDict.items():
        for chrom_section in sorted(chrom_sections):
            of.write('{}\t{}\n'.format(effect, chrom_section))
    of.close()


def OutputFile(geneSegDF, outputDir, filename, includeBands=False):
    """Given a df, output the new file"""
    outputcols = ['genes', 'chr', 'start', 'start_gene', 'start_exon', 'end', 'end_gene', 'segment_end_exon',
                  'Num_Probes', 'sample', 'modal_total_cn', 'expected_total_cn', 'rescaled.cn.a1', 'rescaled.cn.a2',
                  'focality_1', 'focality_2', 'called_CNA1', 'called_CNA2', 'fr_below_1', 'fr_above_1',
                  'fr_below_2', 'fr_above_2']
    if includeBands:
        outputcols.append('band')
        outputcols.append('arm')
        OutputArmSummary(geneSegDF, outputDir, filename)

    geneSegDF.to_csv(os.path.join(outputDir, filename), "\t", header=True, index=False, columns=outputcols)


def getBandInfo(cytoBandDict, chrom, start):
    bands = cytoBandDict['chr{}'.format(chrom)]
    index = np.searchsorted([b['start'] for b in bands], start)
    band = bands[index-1]
    return band


def addCytoBandInfo(geneSegDF):
    package_path = os.path.dirname(os.path.realpath(__file__))
    cytoband_filepath = os.path.join(package_path, 'data/cytoBand_hg19.txt')
    cytoband_df = pd.read_csv(cytoband_filepath, sep='\t', header=None)
    chr_index = 0
    start_index = 1
    band_index = 3
    cytoBandDict = defaultdict(list)
    for index, row in cytoband_df.iterrows():
        chrom = row[chr_index]
        start = row[start_index]
        band = row[band_index]
        cytoBandDict[chrom].append({'start': start, 'band': band, 'arm': band[0]})

    band_info = [getBandInfo(cytoBandDict, row['chr'], row['start']) for index, row in geneSegDF.iterrows()]
    geneSegDF['band'] = pd.Series([b['band'] for b in band_info], index=geneSegDF.index)
    geneSegDF['arm'] = pd.Series([b['arm'] for b in band_info], index=geneSegDF.index)


def main():
    """RUNNING THE CODE"""
    parser = argparse.ArgumentParser(description='Call CNA from ABSOLUTE gene-level segmentation files')
    parser.add_argument('input_dir', metavar='input_dir', type=str)
    parser.add_argument('--build_for_bands', metavar='build_for_bands', type=str)
    parser.add_argument('--output_dir', metavar='output_dir', type=str)
    args = parser.parse_args()

    gene_seg_files_dir = args.input_dir
    output_dir = args.output_dir
    build_for_bands = args.build_for_bands
    if output_dir is None:
        # If no output directory is provided, simply default to placing the output files in the input directory
        output_dir = gene_seg_files_dir

    # Step 1 - grab the list of files
    sys.stdout.write("Getting files in {}\n".format(gene_seg_files_dir))
    annotated_filenames = [n for n in os.listdir(gene_seg_files_dir) if n.endswith('annotated')]

    num_files = len(annotated_filenames)
    if num_files == 0:
        sys.exit("no files found in {}".format(gene_seg_files_dir))

    sys.stdout.write("{} files found in {}\n".format(num_files , gene_seg_files_dir))

    # Step 2 - Iterate through the samples.  For each sample, (1) generate the distribution of revised copy numbers
    # and (2) generate the call (amp, LOH, homozygous deletion, high-level amp).

    # calculate focality = fraction of genome that has rCN > (for deletions) and < (for amplifications)
    # define deletion = rCN < 0.25 and focality > 0.995
    # define amp = focality > 0.98 - 0.2 x log2(rCN/5)
    # define high level amp = focality > 0.98 - (1/7) * log2 (rCN/7)
    count = 1
    for filename in annotated_filenames:  # iterate through each file
        sys.stdout.write("Processing file {} of {}: {}\n".format(str(count), num_files, filename))
        df = pd.read_csv(os.path.join(gene_seg_files_dir, filename), index_col=None, sep='\t', comment='#', low_memory=False)
        CNADist = GenCNADistribution(df)
        GenFocality(df, CNADist)

        # Add hg19 cytobands information column
        if build_for_bands == 'hg19':
            addCytoBandInfo(df)

        OutputFile(df, output_dir, '{}.cna_processed.tsv'.format(filename), includeBands=build_for_bands == 'hg19')
        count += 1

    sys.stdout.write("Processed {} files\n".format(count))


if __name__ == '__main__':
    main()
