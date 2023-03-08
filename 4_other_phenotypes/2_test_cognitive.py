import pandas as pd
import numpy as np

import scipy.stats as stats

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

# Compare the cognitive phenotypes
df=pd.read_csv('analysis_files/1_kendall_cognitive.csv')

# The 2 trail making assessments have low sample size in cases (24), skip these
df.drop(['traila_time', 'trailb_time'], axis=1, inplace=True)

phenos=['pairs_matching_errors','reaction_time_mean', 'fluid_intelligence', 'numeric_memory_maxdigits', 'sd_sub_correctmatches']

# For all others, check skew and normalize, if needed
skew_dat=[]
pdf=PdfPages('Figures/2_normalization.pdf')
for p in phenos:
	arr=df[p].to_list()
	sk=stats.skew(arr, nan_policy='omit', bias=False)
	kurt=stats.kurtosis(arr, nan_policy='omit', bias=False)

	norm_needed=False
	sk2='no normalization'
	kurt2='no normalization'

	if p=='pairs_matching_errors':
		# log+1 normalize
		df[p+'+norm']=df[p].apply(lambda v: np.log10(v+1))
		arr2=df[p+'+norm']
		sk2=stats.skew(arr2, nan_policy='omit', bias=False)
		kurt2=stats.kurtosis(arr2, nan_policy='omit', bias=False)

		norm_needed=True

	if p=='reaction_time_mean':
		# log normalize
		df[p+'+norm']=df[p].apply(lambda v: np.log10(v))
		arr2=df[p+'+norm']
		sk2=stats.skew(arr2, nan_policy='omit', bias=False)
		kurt2=stats.kurtosis(arr2, nan_policy='omit', bias=False)

		norm_needed=True

	if p=='sd_sub_correctmatches':
		# Remove outliers (<3 and >36 subsitutions)
		df[p+'+norm']=df[p]
		df.loc[(df[p+'+norm']<3) | (df[p+'+norm']>36), p+'+norm']=np.nan
		arr2=df[p+'+norm']
		sk2=stats.skew(arr2, nan_policy='omit', bias=False)
		kurt2=stats.kurtosis(arr2, nan_policy='omit', bias=False)

		norm_needed=True

	# Plot before and after distributions
	fig, ax1=plt.subplots()
	sns.histplot(arr, color='blue', alpha=0.5, ax=ax1)
	ax1.tick_params(axis='x', color='blue')

	if norm_needed:
		ax2=ax1.twiny()
		sns.histplot(arr2, color='red', alpha=0.5, ax=ax2)
		ax2.tick_params(axis='x', color='red')

	plt.title

	pdf.savefig()
	plt.close()

	skew_dat.append([p, sk, sk2, kurt, kurt2])
pdf.close()
skew_df=pd.DataFrame(skew_dat, columns=['Phenotype', 'Skew', 'Skew_after_norm', 'Kurtosis', 'Kurtosis_after_norm'])
# Save skew data
skew_df.to_csv('result_tables/2_cognition_skew.csv', index=False)


# Convert to Z score
score_cols=['pairs_matching_errors+norm','reaction_time_mean+norm', 'fluid_intelligence', 'numeric_memory_maxdigits', 'sd_sub_correctmatches+norm']
for sc in score_cols:
	name=sc.split('+')[0]

	vals=df[~df[sc].isnull()][sc].to_numpy()

	score_mean=np.mean(vals)
	score_sd=np.std(vals)

	df[name+'_Z']=(df[sc] - score_mean)/score_sd

	# Change direction of scores so that lower score is worse and higher score is better
	if sc in ['pairs_matching_errors+norm','reaction_time_mean+norm']:
		df[name+'_Z']=df[name+'_Z']*-1


# Perform t-tests and plot distributions
z_cols = [i for i in df.columns.to_list() if 'Z' in i]
controls=['NoCNV_Control', 'LargeRare_Control', 'NEJM_Control']
stat_lst=[]
pdf=PdfPages('Figures/2_cognition_histograms.pdf')
for z in z_cols:
	for c in controls:
		case_vals=df[(df.Case_Control=='Case') & (~df[z].isnull())][z].to_numpy()
		cont_vals=df[(df.Case_Control==c) & (~df[z].isnull())][z].to_numpy()

		# T-test
		t, p = stats.ttest_ind(case_vals, cont_vals, alternative='less')
		stat_lst.append([z, c, 'one tailed t-test, less', len(case_vals), len(cont_vals), np.nanmean(case_vals), np.nanmean(cont_vals), t, p])

		# Histograms
		sns.histplot(data=df[df.Case_Control.isin(['Case', c])], x=z, hue='Case_Control', hue_order=['Case', c], alpha=0.5)
		pdf.savefig()
		plt.close()

pdf.close()

stat_df=pd.DataFrame(stat_lst, columns=['Measure', 'Control', 'Test', 'Case_n', 'Control_n', 'Case_mean', 'Control_mean', 'statistic', 'p-value'])
print(stat_df)
stat_df.to_csv('result_tables/2_cognition_ttest.csv', index=False)
