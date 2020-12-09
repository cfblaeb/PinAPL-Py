#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 17:15:19 2019

@author: philipp
"""

# =======================================================================
# Sort genes by accumulated fold-change
# =======================================================================
import pandas
import numpy as np
from statsmodels.distributions.empirical_distribution import ECDF


def compute_SigmaFC(sgRNAList, config):
	# ------------------------------------------------
	# Compute SigmaFC Score (non-parallel)
	# ------------------------------------------------
	print('Computing accumulated fold-change scores ...')
	df = sgRNAList.copy()  # maybe not required....
	df[config['GeneMetric']] = np.log10((df['counts']+config['delta'])/(df['control mean']+config['delta']))
	sgRNA_per_gene = df.groupby('gene')['sgRNA'].count()
	sgRNA_per_gene.name = '# sgRNAs'
	df2 = pandas.merge(df.groupby('gene')[config['GeneMetric']].sum(), sgRNA_per_gene, left_index=True, right_index=True)

	sig_grnas_per_gene = df.groupby('gene')['significant'].sum()
	sig_grnas_per_gene.name = "# signif. sgRNAs"

	df2 = df2.merge(sig_grnas_per_gene, left_index=True, right_index=True)
	df2[config['GeneMetric']] *= df2['# signif. sgRNAs']

	# ------------------------------------------------
	# Estimate SigmaFC Score distribution (Permutation)
	# ------------------------------------------------

	print(f'Estimating null distribution ({config["Np"]} permutations)...')
	P_set = np.random.choice(len(sgRNAList), size=(config['Np'], int(df2['# sgRNAs'].median())), replace=True)
	SigmaFC_null = [df.loc[P][config['GeneMetric']].sum() * df.loc[P].significant.sum() for P in P_set]
	ecdf = ECDF(SigmaFC_null, side='left')
	if config['ScreenType'] == 'enrichment':
		sigs = df2[df2['# sgRNAs'] > 1][config['GeneMetric']].apply(lambda x: 1 - ecdf(x))
	else:  # ScreenType == 'depletion':
		sigs = df2[df2['# sgRNAs'] > 1][config['GeneMetric']].apply(lambda x: ecdf(x))
	sigs.name = "p_value"
	df2 = df2.merge(sigs, how='left', left_index=True, right_index=True)
	df2['significant'] = df2['p_value'].apply(lambda x: x < config['alpha_g'])

	return df2.reset_index()
