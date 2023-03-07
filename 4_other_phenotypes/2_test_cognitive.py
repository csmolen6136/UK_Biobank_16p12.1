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

# For all others, check skew and log+1 normalize
skew_dat=[]
pdf=PdfPages('Figures/2_normalization.pdf')
controls=['NoCNV_Control', 'LargeRare_Control', 'NEJM_Control']
for p in df.columns.to_list():
	if p=='Sample' or p=='Case_Control':
		continue
	arr=df[p].to_list()
	sk=stats.skew(arr, nan_policy='omit', bias=False)
	kurt=stats.kurtosis(arr, nan_policy='omit', bias=False)

	df[p+'+norm']=df[p].apply(lambda v: np.log10(v+1))
	arr2=df[p+'+norm']
	sk2=stats.skew(arr2, nan_policy='omit', bias=False)
	kurt2=stats.kurtosis(arr2, nan_policy='omit', bias=False)

	# Plot before and after distributions
	fig, ax1=plt.subplots()
	sns.histplot(arr, color='blue', alpha=0.5, ax=ax1)
	ax1.tick_params(axis='x', color='blue')

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
