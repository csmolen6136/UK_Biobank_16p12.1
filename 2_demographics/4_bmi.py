import pandas as pd

from scipy import stats

# Look for differences in age at death between cases and controls
df=pd.read_csv('analysis_files/2_additional_phenotypes.csv')
df=df[~df['BMI'].isnull()][['eid', 'BMI', 'Case_Control']]
print(df)

# T-tests
stat_lst=[]
case_vals=df[df.Case_Control=='Case']['BMI'].to_numpy()
controls=['NoCNV_Control', 'LargeRare_Control', 'NEJM_Control']
for c in controls:
	control_vals=df[df.Case_Control==c]['BMI'].to_numpy()
	stat, p = stats.ttest_ind(case_vals, control_vals, alternative='greater')
	stat_lst.append([c, 'one tailed t-test, greater', len(case_vals), len(control_vals), stat, p])

# Save stats to file
stat_df=pd.DataFrame(stat_lst, columns=['control', 'tets', 'case_n', 'control_n', 'statistic', 'p'])
print(stat_df)
stat_df.to_csv('result_tables/4_bmi.csv', index=False)

# Plot
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

pdf=PdfPages('figures/4_bmi.pdf')
for c in controls:
	subdf=df[df.Case_Control.isin([c, 'Case'])]
	sns.kdeplot(data=subdf, x='BMI', hue='Case_Control', hue_order=['Case', c], common_norm=False)
	plt.title('BMI - '+c)
	pdf.savefig()
	plt.close()
pdf.close()
