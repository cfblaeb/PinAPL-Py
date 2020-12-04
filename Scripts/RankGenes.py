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
	#r = config['NumGuidesPerGene']

	# ------------------------------------------------
	# Read sgRNA enrichment/depletion table
	# ------------------------------------------------
	os.chdir(sgRNARanksDir)
	print('Loading sgRNA ' + screentype + ' table ...')
	filename = glob(sample + '_*sgRNAList.txt')[0]
	sgRNARanking = pandas.read_csv(filename, sep='\t')
	#if screentype == 'enrichment':
	#	sgRNARanking = sgRNARanking.sort_values(['significant', 'p-value', 'fold change', 'sgRNA'], ascending=[0, 1, 0, 1])
	#elif screentype == 'depletion':
	#	sgRNARanking = sgRNARanking.sort_values(['significant', 'p-value', 'fold change', 'sgRNA'], ascending=[0, 1, 1, 1])
	#genes = list(sgRNARanking['gene'])
	#geneList = list(set(genes))

	#NB_pval = list(sgRNARanking['p-value'])
	#G = len(geneList)
	#NB_sig = list(sgRNARanking['significant'])

	# ------------------------------------------------
	# Find number of sgRNAs per gene (some genes have less than r sgRNAs)
	# ------------------------------------------------
	#nGuides = sgRNARanking.groupby('gene')['sgRNA'].count()

	# ------------------------------------------
	# Find number of significant sgRNAs per gene
	# ------------------------------------------
	# Find genes with at least 1 signif. sgRNA
	print('Looking for sgRNAs with significant fold change ...')
	#Temp_DF = pandas.DataFrame(data={'gene': genes, 'NB_significant': NB_sig}, columns=['gene', 'NB_significant'])
	#Temp_DF = Temp_DF.sort_values(['NB_significant'], ascending=False)
	#NB_sig = list(Temp_DF['NB_significant'])
	#genes0 = list(Temp_DF['gene'])
	#n0 = NB_sig.index(False)
	#sigGenes = [genes0[k] for k in range(n0)]  # genes with at least 1 signif. sgRNA
	#sigGenes = sgRNARanking[sgRNARanking.significant].gene.unique()
	#sigGenesList = list(set(sigGenes))
	# count significant sgRNAs for all genes with at least 1 sign. sgRNA
	#guidesPerGene = list()
	#GeneCounts = Counter()
	#for gene in sigGenes:
	#	GeneCounts[gene] += 1
	#for gene in sigGenesList:
	#	guidesPerGene.append(GeneCounts[gene])

	# Save number of significant sgRNAs per gene
	#sigGuides = list()
	#for gene in geneList:
	#	if gene not in sigGenesList:
	#		sigGuides.append(0)
	#	else:
	#		sigGuides.append(guidesPerGene[sigGenesList.index(gene)])

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
	#os.chdir(EffDir)

	#if len(guidesPerGene) > 0:
	#	plt.hist(guidesPerGene, bins=range(1, r + 2), align='left', color='#42f4a1', edgecolor='black')
	#else:
	#	plt.figtext(0.5, 0.5, 'N/A')
	#plt.title('sgRNA on-Target Efficacy', fontsize=12)
	#plt.xlabel('# Significant sgRNAs', fontsize=11)
	#plt.ylabel('Number of Genes', fontsize=11)
	#plt.xticks(range(1, r + 1))
	#plt.tick_params(labelsize=11)
	#plt.tight_layout()
	#plt.savefig(sample + '_sgRNA_Efficacy.png', dpi=res)
	#if svg:
	#	plt.savefig(sample + '_sgRNA_Efficacy.svg')

	# ------------------------------------------
	# Rank genes according to specified method
	# ------------------------------------------
	os.chdir(ScriptsDir)
	if GeneMetric == 'SigmaFC':
		metric, metric_pval, metric_sig = compute_SigmaFC(sgRNARanking)
		SortFlag = False if screentype == 'enrichment' else True  # metric based on fold-change
	elif GeneMetric == 'aRRA':
		metric, metric_pval, metric_sig = compute_aRRA(sgRNARanking, config)
		SortFlag = True  # metric based on p-val
	elif GeneMetric == 'STARS':
		metric, metric_pval, metric_sig = compute_STARS(sgRNARanking)
		SortFlag = False  # metric always from high to low
	elif GeneMetric == 'AvgLogFC':
		Results_df = compute_avg_log_fc(sgRNARanking, config)
		SortFlag = False if screentype == 'enrichment' else True  # metric based on fold-change
	else:
		print('### ERROR: Cannot find gene ranking method! ###')
	# Correcting cdf artifact in case of no significant sgRNAs
	if Results_df.significant.all():
		Results_df.significant = Results_df.significant == False

	# -------------------------------------------------
	# Plotting p-value distribution
	# -------------------------------------------------
	if not Results_df.p_value.isna().all():
		print('Plotting gene metric p-value distribution...')
		PlotTitle = 'Gene ' + screentype.capitalize() + ' (' + GeneMetric + ')'
		pvalHist(Results_df.p_value, pvalDir, sample, res, svg, '#c8d1ca', PlotTitle)

	# -------------------------------------------------
	# Output list
	# -------------------------------------------------
	if not os.path.exists(GeneDir):
		os.makedirs(GeneDir)
	os.chdir(GeneDir)
	print('Writing results dataframe ...')
	#Results_df = pandas.DataFrame(data={'gene': [geneList[g] for g in range(G)],
	#									GeneMetric: [metric[g] for g in range(G)],
	#									'p_value': [metric_pval[g] for g in range(G)],
	#									'significant': [str(metric_sig[g]) for g in range(G)],
	#									'# sgRNAs': [nGuides[g] for g in range(G)],
	#									'# signif. sgRNAs': [sigGuides[g] for g in range(G)]},
	#							  columns=['gene', GeneMetric, 'p_value', 'significant', '# sgRNAs', '# signif. sgRNAs'])
	Results_df_0 = Results_df.sort_values(['significant', GeneMetric], ascending=[False, SortFlag])
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
