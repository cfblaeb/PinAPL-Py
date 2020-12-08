#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 17:15:19 2019

@author: philipp
"""

# =======================================================================
# Sort genes by adjusted robust rank aggregation (Li et al., Genome Biology 2014)
# =======================================================================
from scipy.stats import beta
from statsmodels.distributions.empirical_distribution import ECDF
from pandas import DataFrame
from numpy.random import choice


def aRRA(g, geneList, genes_x, NB_pval_x, P0_aRRA, ranks_x, L):
	gene = geneList[g]
	gene_ranks_sig = list()
	I = genes_x.index(gene)
	i = I
	terminate = False
	while genes_x[i] == gene and terminate == False:
		if NB_pval_x[i] < P0_aRRA:
			gene_ranks_sig.append(ranks_x[i])
		if i <= L - 2:
			i += 1
		else:
			terminate = True
	gene_ranks_sig.sort()
	J = len(gene_ranks_sig)
	if J > 0:
		U = [gene_ranks_sig[i] / L for i in range(J)]
		rho_k = list()
		for k in range(J):
			kth_smallest = U[k]
			pval = beta.cdf(kth_smallest, k + 1, J - k)
			rho_k.append(pval)
		rho = min(rho_k)
	else:  # no sgRNAs are significant
		rho = 1
	return rho


def aRRA_null(I, NB_pval_x, P0_aRRA, ranks, L):
	I_ranks_sig = list()
	for i in I:
		if NB_pval_x[i] < P0_aRRA:
			I_ranks_sig.append(ranks[i])
	I_ranks_sig.sort()
	J = len(I_ranks_sig)
	if J > 0:
		U = [I_ranks_sig[i] / L for i in range(J)]
		rho_k = list()
		for k in range(J):
			kth_smallest = U[k]
			pval = beta.cdf(kth_smallest, k + 1, J - k)
			rho_k.append(pval)
		rho_null = min(rho_k)
	else:
		rho_null = 1
	return rho_null


def compute_aRRA(HitList, config):
	Np = config['Np']
	r = config['NumGuidesPerGene']
	alpha_g = config['alpha_g']
	screentype = config['ScreenType']
	fc = list(HitList['fold change'])
	NB_pval = list(HitList['p-value'])
	sig = list(HitList['significant'])
	genes = list(HitList['gene'])
	L = len(genes)
	geneList = list(set(genes))
	G = len(geneList)

	# -------------------------------------------
	# Compute fold change ranks
	# -------------------------------------------
	fc_DF = DataFrame(data={'fc': fc}, columns=['fc'])
	if screentype == 'enrichment':
		fc_DF = fc_DF.rank(ascending=False)
	elif screentype == 'depletion':
		fc_DF = fc_DF.rank(ascending=True)

	ranks = list(fc_DF['fc'])
	# -------------------------------------------
	# Read out maximal sgRNA p-value that's still significant
	# -------------------------------------------
	I0 = sig.index(False)
	P0_aRRA = NB_pval[I0 - 1]
	# -------------------------------------------
	# Compute aRRA metric
	# -------------------------------------------
	if set(['N/A']) != set(NB_pval):
		print('Computing aRRA scores ...')
		aRRA_DF = DataFrame(data={'gene': [genes[i] for i in range(L)], 'ranks': [ranks[i] for i in range(L)], 'NB_pval': [NB_pval[i] for i in range(L)]}, columns=['gene', 'ranks', 'NB_pval'])
		aRRA_DF = aRRA_DF.sort_values(['gene'])
		genes_x = list(aRRA_DF['gene'])
		ranks_x = list(aRRA_DF['ranks'])
		NB_pval_x = list(aRRA_DF['NB_pval'])
		metric = [aRRA(g, geneList, genes_x, NB_pval_x, P0_aRRA, ranks_x, L) for g in range(G)]
		# Permutation
		print('Estimating aRRA null distribution (' + str(Np) + ' permutations)...')
		I_perm = choice(L, size=(Np, r), replace=True)
		metric_null = [aRRA_null(I, NB_pval_x, P0_aRRA, ranks, L) for I in I_perm]
		# p-value
		print('Estimating aRRA p-values ...')
		ecdf = ECDF(metric_null)
		metric_pval = list()
		for g in range(G):
			pval = ecdf(metric[g])
			metric_pval.append(pval)
		sig_metric = [True if metric_pval[g] < alpha_g else False for g in range(G)]
	else:  # no control replicates
		print('### ERROR: Cannot compute aRRA scores without significant sgRNAs! ###')
		metric = [-1 for k in range(G)]
		metric_pval = [-1 for k in range(G)]

	sgRNA_series = HitList.groupby('gene')['sgRNA'].count()
	sgRNA_series.name = '# sgRNAs'

	sig_grnas_per_gene = HitList.groupby('gene')['significant'].sum()
	sig_grnas_per_gene.name = "# signif. sgRNAs"

	df = DataFrame({'gene': geneList, config['GeneMetric']: metric, 'p_value': metric_pval, 'significant': sig_metric})
	df = df.merge(sgRNA_series, left_on='gene', right_index=True)
	df = df.merge(sig_grnas_per_gene, left_on='gene', right_index=True)
	return df
