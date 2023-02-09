import pandas as pd

from scipy import stats

# Look for differences in age at death between cases and controls
df=pd.read_csv('analysis_files/2_additional_phenotypes.csv')
imd_cols=[i for i in df.columns.to_list() if 'England' in i]
df.dropna(subset=imd_cols, how='all', inplace=True)
df=df[['eid', 'Case_Control']+imd_cols]
print(df)

# T-tests
stat_lst=[]
controls=['NoCNV_Control', 'LargeRare_Control', 'NEJM_Control']
for imd in imd_cols:
	case_vals=df[df.Case_Control=='Case'][imd].to_numpy()
	for c in controls:
		control_vals=df[df.Case_Control==c][imd].to_numpy()
		stat, p = stats.ttest_ind(case_vals, control_vals, alternative='greater')
		stat_lst.append([imd, c, 'one tailed t-test, greater', len(case_vals), len(control_vals), stat, p])

# Save stats to file
stat_df=pd.DataFrame(stat_lst, columns=['measure', 'control', 'test', 'case_n', 'control_n', 'statistic', 'p'])
print(stat_df)
stat_df.to_csv('result_tables/5_imd.csv', index=False)

# Plot
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

pdf=PdfPages('figures/5_imd.pdf')
for imd in imd_cols:
	for c in controls:
		subdf=df[df.Case_Control.isin([c, 'Case'])]
		sns.kdeplot(data=subdf, x=imd, hue='Case_Control', hue_order=['Case', c], common_norm=False)
		plt.title(imd+' - '+c)
		pdf.savefig()
		plt.close()
pdf.close()
