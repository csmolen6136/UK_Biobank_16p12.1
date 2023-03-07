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
print(skew_df)
