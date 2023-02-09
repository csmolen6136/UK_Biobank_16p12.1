import pandas as pd

# Clean up the results from the logistic models
df=pd.read_csv('result_tables/5_logistic_models.csv')

# In R, the node IDs had Xs added to the front - remove
df=df[df.id!='.']
df['id']=df['id'].str.split('X', expand=True)[1].astype(int)

# Bonferroni correct by control and category
controls=['NoCNV_Control', 'LargeRare_Control', 'NEJM_Control']
types=['node', 'Sub-block', 'Block', 'Chapter']

count_dict={}
for c in controls:
	for t in types:
		count=df[(df.type==t) & (df.control==c)].shape[0]
		# Multiply by 2 - case/control and sex
		count_dict[c+'.'+t]=count*2
df['cont_type']=df.control+'.'+df.type
df['factor']=df.cont_type.map(count_dict)

df['cc_bonferroni_fdr']=df.cc_p*df.factor
df.loc[df.cc_bonferroni_fdr>1, 'cc_bonferroni_fdr']=1
df['sex_bonferroni_fdr']=df.sex_p*df.factor
df.loc[df.sex_bonferroni_fdr>1, 'sex_bonferroni_fdr']=1

# Add pos/neg to indicate direction
df['cc_posneg']='-'
df.loc[df.cc_est>=0, 'cc_posneg']='+'
df['sex_posneg']='-'
df.loc[df.sex_est>=0, 'sex_posneg']='+'

# Add node meaning
coding=pd.read_csv('analysis_files/coding19_updated.csv')
df['meaning']=df['id'].map(dict(zip(coding.node_id.to_list(), coding.meaning.to_list())))

# Remove extra columns
df=df[['type', 'id', 'control', 'cc_est', 'cc_err', 'cc_p', 'cc_bonferroni_fdr', 'cc_posneg',
	'sex_est', 'sex_err', 'sex_p', 'sex_bonferroni_fdr', 'sex_posneg',
	'meaning']]

print(df)

for c in controls:
	print(df[(df.control==c) & (df.cc_bonferroni_fdr<=0.05)])


# Save to file
df.to_csv('result_tables/6_parsed_model.csv', index=False)
