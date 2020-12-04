#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 17:15:19 2019

@author: philipp
"""

# =======================================================================
# Sort genes by average Log FC 
# =======================================================================
import numpy
from statsmodels.distributions.empirical_distribution import ECDF


#def AverageLogFC(g, geneList, GeneBundle, lfc, L, repl_avg):
#	gene = geneList[g]
#	lfc_list = list()
#	i0 = GeneBundle.index(gene)
#	i = i0
#	terminate = False
#	while GeneBundle[i] == gene and terminate == False:
#		lfc_list.append(lfc[i])
#		if i <= L - 2:
#			i += 1
#		else:
#			terminate = True
#	if repl_avg == 'mean':
#		avglfc = numpy.mean(lfc_list)
#	elif repl_avg == 'median':
#		avglfc = numpy.median(lfc_list)
#	return avglfc


def AvgLogFC_null(I, lfc, repl_avg):
	logFC_I = [lfc[i] for i in I]
	if repl_avg == 'mean':
		avglfc_I = numpy.mean(logFC_I)
	else:
		avglfc_I = numpy.median(logFC_I)
	return avglfc_I


def compute_AvgLogFC(HitList, nGuides, config):
	r = config['NumGuidesPerGene']
	screentype = config['ScreenType']
	Np = config['Np']
	alpha_g = config['alpha_g']
	repl_avg = config['repl_avg']
	# -------------------------------------------------
	# Compute average log fold change across sgRNAs
	# -------------------------------------------------
	print('Computing average fold change across sgRNAs ...')
	#Aux_DF = pandas.DataFrame(data={'gene': [genes[i] for i in range(L)], 'lfc': [numpy.log2(fc[i]) for i in range(L)]}, columns=['gene', 'lfc'])
	#Aux_DF = Aux_DF.sort_values(['gene'])
	HitList['lfc'] = HitList['fold change'].apply(numpy.log2)
	HitList.sort_values(by='gene', inplace=True)

	lfc = list(HitList['lfc'])
	if repl_avg == 'mean':
		AvgLogFC = HitList.groupby('gene')['lfc'].mean()[set(HitList.gene)].values
	else:
		AvgLogFC = HitList.groupby('gene')['lfc'].median()[set(HitList.gene)].values
	#AvgLogFC = [AverageLogFC(g, geneList, GeneBundle, lfc, L, repl_avg) for g in range(G)]
	# -------------------------------------------------
	# Compute permutations
	# -------------------------------------------------
	I_perm = numpy.random.choice(len(HitList), size=(Np, r), replace=True)
	metric_null = [AvgLogFC_null(I, lfc, repl_avg) for I in I_perm]
	ecdf = ECDF(metric_null)
	metric_pval = list()
	for g in range(len(set(HitList['gene']))):
		if nGuides[g] == 1:
			pval = 'N/A'  # exclude p-value for miRNAs etc
			metric_pval.append(pval)
		elif screentype == 'enrichment':
			pval = 1 - ecdf(AvgLogFC[g])
			metric_pval.append(pval)
		elif screentype == 'depletion':
			pval = ecdf(AvgLogFC[g])
			metric_pval.append(pval)
		else:
			print('### ERROR: Check spelling of ScreenType in configuration file! ###')
	metric_sig = [True if isinstance(metric_pval[g], float) and metric_pval[g] < alpha_g else False for g in range(len(set(HitList['gene'])))]
	return AvgLogFC, metric_pval, metric_sig
