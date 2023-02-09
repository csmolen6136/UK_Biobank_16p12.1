import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Plot the results of the logistic models
df=pd.read_csv('result_tables/6_parsed_model.csv')

# Only plot significant results
df=df[df.cc_bonferroni_fdr<=0.05]

# Plot results by control, by type
pdf=PdfPages('figures/7_logistic_model_results.pdf')
controls=['NoCNV_Control', 'LargeRare_Control', 'NEJM_Control']
types=['node', 'Sub-block', 'Block', 'Chapter']
for c in controls:
	for t in types:
		subdf=df[(df.control==c) & (df.type==t)]
		if subdf.shape[0]==0:
			continue
		# Make a forest plot showing all significant results
		ys=list(range(0, subdf.shape[0]))
		# Plot Case-Control estimates and error bars
		plt.scatter(subdf.cc_est.to_list(), [i+0.1 for i in ys], color='#2c3e50', label='Case-Control')
		for i in ys:
			est=subdf.cc_est.to_list()[i]
			err=subdf.cc_err.to_list()[i]
			plt.plot([est-err, est+err], [i+0.1, i+0.1], color='#2c3e50')
		# Plot sex estimates and error bars
		plt.scatter(subdf.sex_est.to_list(), [i-0.1 for i in ys], color='#f1c40f', label='Sex')
		for i in ys:
			est=subdf.sex_est.to_list()[i]
			err=subdf.sex_err.to_list()[i]
			plt.plot([est-err, est+err], [i-0.1, i-0.1], color='#f1c40f')
		# Add axis labels
		meanings=subdf.meaning.to_list()
		ylabels=[]
		for i in meanings:
			if len(i)>60:
				ylabels.append(i[0:60]+' . . .')
			else:
				ylabels.append(i)
		plt.yticks(ys, ylabels)
		# Add title
		plt.title(c+' '+t+' significant results')
		# Add a line at 0
		plt.plot([0, 0], [-0.5, max(ys)+0.5], color='k', ls=':', zorder=0)
		plt.ylim(-0.25, max(ys)+0.25)
		# Add a legend
		plt.legend()

		plt.tight_layout()
		pdf.savefig()
		plt.close()
pdf.close()
print(df)
