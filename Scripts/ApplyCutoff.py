#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 15:54:46 2018

@author: philipp
"""

# Apply (optional) read count cutoff and obtain counts per gene
# =======================================================================
# Imports
import os
import pandas


def ReadCountCutoff(sample, config):
	# ------------------------------------------------
	# Get parameters
	# ------------------------------------------------
	sgRNAReadCountDir = config['sgRNAReadCountDir']
	GeneReadCountDir = config['GeneReadCountDir']
	minN = config['Cutoff']
	GuideCount_Suffix = '_GuideCounts.txt'
	GeneCount_Suffix = '_GeneCounts.txt'

	# ------------------------------------------------
	# Read library
	# ------------------------------------------------

	# Load sgRNA counts
	os.chdir(sgRNAReadCountDir)
	ReadCountTable = pandas.read_csv(sample + '_GuideCounts.txt', sep='\t', names=['ID', 'gene', 'counts'])

	# Apply read count cut-off
	if minN > 0:
		print('Applying sgRNA minimal count cutoff...')

	# Write output
	# sgRNAs
	GuideCountsFilename = sample + GuideCount_Suffix
	ReadCountTable.to_csv(GuideCountsFilename, sep='\t', index=False, header=False)
	# genes
	if not os.path.exists(GeneReadCountDir):
		os.makedirs(GeneReadCountDir)
	os.chdir(GeneReadCountDir)
	GeneCountsFilename = sample + GeneCount_Suffix
	GeneCounts = ReadCountTable.groupby('gene')['counts'].sum().reset_index().sort_values('gene')
	GeneCounts.to_csv(GeneCountsFilename, sep='\t', index=False, header=False)
