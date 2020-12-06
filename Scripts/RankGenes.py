#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 14:46:37 2016

@author: philipp
"""

# Find candidate genes
# =======================================================================
# Imports
import os
import time
import pandas
from glob import glob
from matplotlib import pyplot as plt
from pvalPlots import pvalHist
from RankGenes_SigmaFC import compute_SigmaFC
from RankGenes_AvgLogFC import compute_avg_log_fc
from RankGenes_aRRA import compute_aRRA
from RankGenes_STARS import compute_STARS


def GeneRankingAnalysis(sample, config):
	# ------------------------------------------------
	# Print header
	# ------------------------------------------------
	print('++++++++++++++++++++++++++++++++++++++++++++++++')
	start_total = time.time()

	# ------------------------------------------------
	# Get parameters
	# ------------------------------------------------
	ScriptsDir = config['ScriptsDir']
	sgRNARanksDir = config['sgRNARanksDir']
	EffDir = config['EffDir']
	GeneDir = config['GeneDir']
	screentype = config['ScreenType']
	GeneMetric = config['GeneMetric']
	SheetFormat = config['HitListFormat']
	pvalDir = config['pvalDir_genes']
	res = config['dpi']
	svg = config['svg']

	# ------------------------------------------------
	# Read sgRNA enrichment/depletion table
	# ------------------------------------------------
	os.chdir(sgRNARanksDir)
	print('Loading sgRNA ' + screentype + ' table ...')
	filename = glob(sample + '_*sgRNAList.txt')[0]
	sgRNARanking = pandas.read_csv(filename, sep='\t')

	# Plot histogram
	if not os.path.exists(EffDir):
		os.makedirs(EffDir)
	f, ax = plt.subplots(figsize=(3.5, 2.9))
	sgRNARanking.groupby('gene')['significant'].sum().value_counts().sort_index().plot(kind='bar', ax=ax)
	ax.set_title('sgRNA on-Target Efficacy', fontsize=12)
	ax.set_xlabel('# Significant sgRNAs', fontsize=11)
	ax.set_ylabel('Number of Genes', fontsize=11)
	plt.savefig(EffDir + "/" + sample + '_sgRNA_Efficacy.png', dpi=res)
	if svg:
		f.savefig(EffDir + "/" + sample + '_sgRNA_Efficacy.svg')

	# ------------------------------------------
	# Rank genes according to specified method
	# ------------------------------------------
	os.chdir(ScriptsDir)
	if GeneMetric == 'SigmaFC':
		results_df = compute_SigmaFC(sgRNARanking, config)
		sort_flag = False if screentype == 'enrichment' else True  # metric based on fold-change
	elif GeneMetric == 'aRRA':
		results_df = compute_aRRA(sgRNARanking, config)
		sort_flag = True  # metric based on p-val
	elif GeneMetric == 'STARS':
		results_df = compute_STARS(sgRNARanking, config)
		sort_flag = False  # metric always from high to low
	else:  # GeneMetric == 'AvgLogFC':
		results_df = compute_avg_log_fc(sgRNARanking, config)
		sort_flag = False if screentype == 'enrichment' else True  # metric based on fold-change

	# Correcting cdf artifact in case of no significant sgRNAs
	if results_df.significant.all():
		results_df.significant = results_df.significant == False

	# -------------------------------------------------
	# Plotting p-value distribution
	# -------------------------------------------------
	if not results_df.p_value.isna().all():
		print('Plotting gene metric p-value distribution...')
		PlotTitle = 'Gene ' + screentype.capitalize() + ' (' + GeneMetric + ')'
		pvalHist(results_df.p_value, pvalDir, sample, res, svg, '#c8d1ca', PlotTitle)

	# -------------------------------------------------
	# Output list
	# -------------------------------------------------
	if not os.path.exists(GeneDir):
		os.makedirs(GeneDir)
	os.chdir(GeneDir)
	print('Writing results dataframe ...')
	Results_df_0 = results_df.sort_values(['significant', GeneMetric], ascending=[False, sort_flag])
	GeneListFilename = filename[0:-14] + '_' + GeneMetric + '_GeneList.txt'
	Results_df_0.to_csv(GeneListFilename, sep='\t', index=False)
	if SheetFormat == 'xlsx':
		print('Converting to xlsx ...')
		GeneListFilename = filename[0:-14] + '_' + GeneMetric + '_GeneList.xlsx'
		Results_df_0.to_excel(GeneListFilename)

	# -------------------------------------------------
	# Time Stamp
	# -------------------------------------------------
	os.chdir(ScriptsDir)
	end_total = time.time()
	print('------------------------------------------------')
	print('Script completed.')
	sec_elapsed = end_total - start_total
	if sec_elapsed < 60:
		time_elapsed = sec_elapsed
		print('Time elapsed [secs]: ' + '%.3f' % time_elapsed + '\n')
	elif sec_elapsed < 3600:
		time_elapsed = sec_elapsed / 60
		print('Time elapsed [mins]: ' + '%.3f' % time_elapsed + '\n')
	else:
		time_elapsed = sec_elapsed / 3600
		print('Time elapsed [hours]: ' + '%.3f' % time_elapsed + '\n')
