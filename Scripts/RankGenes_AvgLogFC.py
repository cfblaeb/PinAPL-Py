"""
Created on Wed May 15 17:15:19 2019

@author: philipp
"""

# =======================================================================
# Sort genes by average Log FC 
# =======================================================================
import numpy
import pandas
from statsmodels.distributions.empirical_distribution import ECDF


def avg_log_fc_null(ind, lfc, repl_avg):
	log_fc_i = [lfc[i] for i in ind]
	if repl_avg == 'mean':
		avg_lfc_i = numpy.mean(log_fc_i)
	else:
		avg_lfc_i = numpy.median(log_fc_i)
	return avg_lfc_i


def compute_avg_log_fc(sgrna_df, config):
	r = config['NumGuidesPerGene']
	screentype = config['ScreenType']
	config_np = config['Np']
	alpha_g = config['alpha_g']
	repl_avg = config['repl_avg']
	# -------------------------------------------------
	# Compute average log fold change across sgRNAs
	# -------------------------------------------------
	print('Computing average fold change across sgRNAs ...')
	sgrna_df['AvgLogFC'] = sgrna_df['fold change'].apply(numpy.log2)

	if repl_avg == 'mean':
		avg_log_fc = sgrna_df.groupby('gene')['AvgLogFC'].mean()
	else:
		avg_log_fc = sgrna_df.groupby('gene')['AvgLogFC'].median()
	no_sgrna = sgrna_df.groupby('gene')['sgRNA'].count()
	no_sgrna.name = "# sgRNAs"
	avg_log_fc = pandas.merge(avg_log_fc, no_sgrna, left_index=True, right_index=True)
	# -------------------------------------------------
	# Compute permutations
	# -------------------------------------------------
	I_perm = numpy.random.choice(len(sgrna_df), size=(config_np, r), replace=True)
	metric_null = [avg_log_fc_null(I, list(sgrna_df['AvgLogFC']), repl_avg) for I in I_perm]
	ecdf = ECDF(metric_null)
	if screentype == 'enrichment':
		sigs = avg_log_fc[avg_log_fc['# sgRNAs'] > 1]['AvgLogFC'].apply(lambda x: 1-ecdf(x))
	elif screentype == 'depletion':
		sigs = avg_log_fc[avg_log_fc['# sgRNAs'] > 1]['AvgLogFC'].apply(lambda x: ecdf(x))
	else:
		raise ValueError(f"screentype not recognized: {screentype}")
	sigs.name = "p_value"
	avg_log_fc = avg_log_fc.merge(sigs, how='left', left_index=True, right_index=True)
	avg_log_fc['significant'] = avg_log_fc['p_value'].apply(lambda x: x < alpha_g)
	sig_grnas_per_gene = sgrna_df.groupby('gene')['significant'].sum()
	sig_grnas_per_gene.name = "# signif. sgRNAs"
	avg_log_fc = avg_log_fc.merge(sig_grnas_per_gene, left_index=True, right_index=True)
	return avg_log_fc
