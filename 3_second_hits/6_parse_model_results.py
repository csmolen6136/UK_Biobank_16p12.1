import pandas as pd
import statsmodels.stats.multitest as stats

# Clean up the results from the logistic models
df=pd.read_csv('result_tables/5_logistic_models.csv')

# Bonferroni correct by control
controls=['NoCNV_Control', 'LargeRare_Control', 'NEJM_Control']

df['cc_bh_fdr']=-1
df['sex_bh_fdr']=-1
for c in controls:
	cc_p=df[df.control==c].cc_p.to_list()
	sex_p=df[df.control==c].sex_p.to_list()
	bh_fdr=list(stats.multipletests(cc_p+sex_p, alpha=0.5, method='fdr_bh')[1])
	print(len(cc_p), len(sex_p), len(bh_fdr))
	df.loc[df.control==c, 'cc_bh_fdr']=bh_fdr[0:len(cc_p)]
	df.loc[df.control==c, 'sex_bh_fdr']=bh_fdr[len(cc_p):]


# Add pos/neg to indicate direction
df['cc_posneg']='-'
df.loc[df.cc_est>=0, 'cc_posneg']='+'
df['sex_posneg']='-'
df.loc[df.sex_est>=0, 'sex_posneg']='+'

# Remove extra columns
df=df[['chrom', 'gene', 'control', 'case_n', 'control_n', 'cc_est', 'cc_err', 'cc_p', 'cc_bh_fdr', 'cc_posneg',
	'sex_est', 'sex_err', 'sex_p', 'sex_bh_fdr', 'sex_posneg']]
print(df)

for c in controls:
	print(df[(df.control==c) & (df.cc_bh_fdr<=0.05)])

# Save to file
df.to_csv('result_tables/6_parsed_model.csv', index=False)
