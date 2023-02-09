import pandas as pd
import numpy as np

import scipy.stats as stats

from statsmodels.stats.proportion import proportions_ztest

# Assess a potential sex bias in 16p12.1 del carriers
df=pd.read_csv('../0_get_carriers/UKB_participant_table.csv')
print(df)

stat_lst=[]
# Difference from 50% - binomial
df_16p=df[df.Case_Control=='Case']
males=sum(df_16p.Sex=='M')

p=stats.binom_test(x=males, n=df_16p.shape[0], p=0.5, alternative='less')
stat_lst.append(['one tailed binomial test - less', '.', males/df_16p.shape[0], 0.5, males/df_16p.shape[0]-0.5, np.nan, p])


# Difference from control groups with Z test
controls=['NoCNV_Control', 'LargeRare_Control', 'NEJM_Control']

case_m=df[(df.Case_Control=='Case') & (df.Sex=='M')].shape[0]
case_tot=df[df.Case_Control=='Case'].shape[0]
for c in controls:
	cont_m=df[(df.Case_Control==c) & (df.Sex=='M')].shape[0]
	cont_tot=df[df.Case_Control==c].shape[0]

	statistic, p = proportions_ztest([case_m, cont_m], [(case_tot), (cont_tot)], alternative='smaller', value=0)

	diff=case_m/case_tot - cont_m/cont_tot

	stat_lst.append(['one tailed Z test - less', c, case_m/case_tot, cont_m/cont_tot, diff, statistic, p])

stat_df=pd.DataFrame(stat_lst, columns=['test', 'control', 'case_proportion', 'control_proportion', 'proportion_difference', 'statistic', 'p'])

# Save to file
stat_df.to_csv('result_tables/1_sex_bias.csv', index=False)

# Plot proportions of males
import matplotlib.pyplot as plt

plt.bar([1, 2, 3, 4], [stat_df.case_proportion.to_list()[0]] + stat_df.control_proportion.to_list()[1:])
plt.plot([0, 5], [0.5, 0.5], color='k', ls=':', alpha=0.5)
plt.xlim(0.5, 4.5)
plt.xticks([1, 2, 3, 4], ['16p12.1 cases', 'No CNV control', 'Large Rare deletion control', 'NEJM del control'])
plt.ylabel('Proportion males')
plt.savefig('figures/1_sex_bias.pdf')
plt.close()
