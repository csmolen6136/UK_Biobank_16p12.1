import pandas as pd
import sys

prefix=int(sys.argv[1])

print('This is run '+str(prefix))
print('It will run terms GO:'+str((prefix-1)*20500).zfill(7)+' - GO:'+str(prefix*20500).zfill(7))

# Create a binary matrix of samples and GO terms
# To save space/time, break dataframe into chunks with a few thousand GO terms per chunk
df=pd.read_csv('analysis_files/11_sample_go_associations.csv')
df=df[['sample', 'go_id']]
df.drop_duplicates(inplace=True)

df['go_num']=df.go_id.str.split(':', expand=True)[1].astype(int)
df=df[(df.go_num<=(prefix*20500)) & (df.go_num>((prefix-1)*20500))]
if df.shape[0]==0:
	print('No terms in this chunk')
	exit()

# Add in sample label
vars=pd.read_csv('analysis_files/2_secondhits.csv')
vars.index=vars.Sample
sex_dict=vars.Sex.to_dict()
cc_dict=vars.Case_Control.to_dict()

go_terms=list(df.go_id.unique())
print(len(go_terms))
outdf=pd.DataFrame(0, index=list(set(vars.Sample.to_list())), columns=go_terms)

# Add in Case/Control status and sex
outdf['Case_Control']=outdf.index.map(cc_dict)
outdf['Sex']=outdf.index.map(sex_dict).map({'F':0, 'M':1})

for g in go_terms:
	idx=go_terms.index(g)
	if idx%1000==0:
		print(idx)
	go_samps=df[df.go_id==g]['sample'].to_list()
	outdf.loc[go_samps, g]=1
print(outdf)

# Save to file
outdf.to_csv('analysis_files/12_go_matrix/'+str(prefix)+'.csv', index=False)
