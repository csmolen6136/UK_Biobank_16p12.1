import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import scipy.stats as stats

# Compare burden between cases and controls
df=pd.read_csv('../0_get_carriers/UKB_participant_table.csv')

pdf=PdfPages('figures/1_burden_kde.pdf')
pdf2=PdfPages('figures/1_burden_violin.pdf')
stat_lst=[]
# Compare burden as kdeplots
hits=['SNV_burden', 'Genes_del', 'Genes_dup']
controls=['NoCNV_Control', 'LargeRare_Control', 'NEJM_Control']
for h in hits:
	for c in controls:
		subdf=df[(~df[h].isnull()) & (df.Case_Control.isin(['Case', c]))]

		# Kdeplot
		sns.kdeplot(data=subdf, x=h, hue='Case_Control', hue_order=['Case', c], cut=0, common_norm=False)
		plt.title(h+' '+c)
		pdf.savefig()
		plt.close()

		# Violinplot
		subdf['split']='.'
		sns.violinplot(data=subdf, y=h, x='split', hue='Case_Control', hue_order=['Case', c], cut=0, split=True, inner='quartiles')
		plt.title(h+' '+c)
		pdf2.savefig()
		plt.close()

		# Compare with t-test
		case_burden=subdf[subdf.Case_Control=='Case'][h].to_numpy()
		control_burden=subdf[subdf.Case_Control!='Case'][h].to_numpy()

		statistic, p = stats.ttest_ind(case_burden, control_burden, alternative='less')
		stat_lst.append(['one tailed t test - less', h, c, len(case_burden), len(control_burden), np.mean(case_burden), np.mean(control_burden), statistic, p])

pdf.close()
pdf2.close()

stat_df=pd.DataFrame(stat_lst, columns=['test', 'trait', 'control', 'case n', 'control n', 'case mean', 'control mean', 'statistic', 'p'])
print(stat_df)

# Save to file
stat_df.to_csv('result_tables/1_burden_ttest.csv', index=False)
