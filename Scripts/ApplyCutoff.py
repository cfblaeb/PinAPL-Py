#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 15:54:46 2018

@author: philipp
"""

# Apply (optional) read count cutoff and obtain counts per gene
# =======================================================================
# Imports
import yaml
import os
import sys
import pandas
#from joblib import Parallel, delayed
#import multiprocessing


#def CountReadsPerGene(g):       # faster than slicing
#    gene = GeneList[g]
#    I = geneIDs.index(gene)
#    j = I
#    g_counts = 0
#    terminate = False
#    while geneIDs[j] == gene and terminate == False:
#        g_counts = g_counts + sgRNA_counts[j]
#        if j <= L-2:
#            j+=1
#        else:
#            terminate = True
#    return g_counts


def ReadCountCutoff(sample):
    # ------------------------------------------------
    # Get parameters
    # ------------------------------------------------
    config = yaml.load(open('configuration.yaml','r'), Loader=yaml.FullLoader)
    sgRNAReadCountDir = config['sgRNAReadCountDir']
    GeneReadCountDir = config['GeneReadCountDir']
    #LibDir = config['LibDir']
    #LibFilename = config['LibFilename']
    #LibFormat = LibFilename[-3:]
    #libsep = '\t' if LibFormat == 'tsv' else ','
    N0 = 1000000
    minN = config['Cutoff']
    GuideCount_Suffix = '_GuideCounts.txt'
    GeneCount_Suffix = '_GeneCounts.txt'    

    # ------------------------------------------------
    # Read library
    # ------------------------------------------------  
    #os.chdir(LibDir)
    #LibCols = ['gene','ID','seq']
    #LibFile = pandas.read_csv(LibFilename, sep=libsep, skiprows=1, names=LibCols)
    #LibFile = LibFile.sort_values(['gene', 'ID'])
    #sgIDs = list(LibFile['ID'])
    #global L
    #L = len(sgIDs)
    #global geneIDs
    #geneIDs = list(LibFile['gene'])
    #G = len(set(geneIDs))

    # Load sgRNA counts
    os.chdir(sgRNAReadCountDir)
    ReadCountTable = pandas.read_csv(sample+'_GuideCounts.txt', sep='\t', names=['ID','gene','counts'])
    ReadsPerGuide = list(ReadCountTable['counts'])
    NReads = ReadCountTable.counts.sum()

    # Apply read count cut-off
    if minN > 0:    
        print('Applying sgRNA minimal count cutoff...')
    ReadSel = [x >= NReads/N0*minN for x in ReadsPerGuide]
    global sgRNA_counts
    sgRNA_counts = [a*b for a,b in zip(ReadSel, ReadsPerGuide)]
    
    # Read counts per gene in library       
    #global GeneList
    #GeneList = list(set(geneIDs))
    #num_cores = multiprocessing.cpu_count()
    #gene_counts = Parallel(n_jobs=num_cores)(delayed(CountReadsPerGene)(g) for g in range(G))
    
    # Write output  
    # sgRNAs
    GuideCountsFilename = sample + GuideCount_Suffix
    #GuideCounts = pandas.DataFrame()
    #GuideCounts['sgIDs'] = sgIDs
    #GuideCounts['geneIDs'] = geneIDs
    #GuideCounts['sgRNA_counts'] = sgRNA_counts
    ReadCountTable.to_csv(GuideCountsFilename, sep = '\t', index = False, header = False)
    # genes
    if not os.path.exists(GeneReadCountDir):
        os.makedirs(GeneReadCountDir)
    os.chdir(GeneReadCountDir)
    GeneCountsFilename = sample + GeneCount_Suffix
    #GeneCounts = pandas.DataFrame()
    #GeneCounts['gene'] = GeneList
    #GeneCounts['gene_counts'] = gene_counts
    #GeneCounts = GeneCounts.sort_values('gene')
    GeneCounts = ReadCountTable.groupby('gene')['counts'].sum().reset_index().sort_values('gene')
    GeneCounts.to_csv(GeneCountsFilename, sep = '\t', index = False, header = False)
    
    
if __name__ == "__main__":
    input1 = sys.argv[1]
    ReadCountCutoff(input1)
