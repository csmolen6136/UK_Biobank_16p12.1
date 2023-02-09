import pandas as pd

# Reformat ICD10 data into table form for a logistic model
# Create a table with each sample as a row and each ICD10 code as a column
df = pd.read_csv('analysis_files/1_icd_meaning.csv')
print(df)

samples=list(df.Sample.unique())

# Create separate files for nodes, Sub-blocks, Blocks, and Chapters
types=['node', 'Sub-block', 'Block', 'Chapter']
# Also add in sample info
table=pd.read_csv('../0_get_carriers/UKB_participant_table.csv')

for type in types:
	nodes=list(df[type].unique())
	node_df=pd.DataFrame(0, index=samples, columns=nodes)
	for n in nodes:
		node_samps=df[df[type]==n]['Sample'].to_list()
		node_df.loc[node_samps, n]=1
	node_df['Sex']=node_df.index.map(dict(zip(table.Sample.to_list(), table.Sex.to_list())))
	node_df.Sex=node_df.Sex.map({'F':0, 'M':1})
	node_df['Case_Control']=node_df.index.map(dict(zip(table.Sample.to_list(), table.Case_Control.to_list())))
	node_df['Sample']=node_df.index.to_list()
	print(node_df)
	# Save to file
	node_df.to_csv('model_files/4_'+type+'.csv', index=False)
