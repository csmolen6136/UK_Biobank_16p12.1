import pandas as pd

# Create a table of samples in the UK Biobank
# Table contains ID, sex, 16p12.1 carrier status, other NEJM CNV carrier status, other large (>500 kb) rare genic CNV carrier status, ICD10 code availability, 2nd hits availability
# This script will only add in CNV data

cnvs=pd.read_csv('/data5/UK_Biobank/PennCNV_calls/bed_files/19_gene_filter.bed', sep='\t')
samples=list(cnvs.Sample.unique())

df=pd.DataFrame({'Sample':samples, 'CNV_data_available':['X']*len(samples)})
print(df)

# Add samples with 16p12.1
samples16p=cnvs[(cnvs.NEJM_Name.str.contains('16p12.1')) & (cnvs.Type=='DEL')]['Sample'].to_list()
df['16p12.1del_carrier']='.'
df.loc[df.Sample.isin(samples16p), '16p12.1del_carrier']='X'
print(df)
print(df['16p12.1del_carrier'].value_counts())

# Add 16p12.1 breakpoints
breakpoints={}
for s in samples16p:
	subcnv=cnvs[(cnvs.Sample==s) & (cnvs.NEJM_Name.str.contains('16p12.1')) & (cnvs.Type=='DEL')]
	points='chr16:'+str(subcnv.Start.to_list()[0])+'-'+str(subcnv.End.to_list()[0])
	breakpoints[s]=points
df['16p12.1del_breakpoints']=df.Sample.map(breakpoints)
df.fillna('.', inplace=True)

# Add other CNVs
othernejm_del=cnvs[(cnvs.NEJM_Name!='.') & ~(cnvs.NEJM_Name.str.contains('16p12.1')) & (cnvs.Type=='DEL')]['Sample'].to_list()
df['NEJM_del_carrier']='.'
df.loc[df.Sample.isin(othernejm_del), 'NEJM_del_carrier']='X'
# Add in specific CNV names
nejmdel_dict={}
for s in othernejm_del:
	names=[i for i in list(cnvs[(cnvs.Sample==s) & (cnvs.Type=='DEL')]['NEJM_Name'].unique()) if i!='.']
	nejmdel_dict[s]=';'.join(names)
df['NEJM_del_name']=df.Sample.map(nejmdel_dict)
df.fillna('.', inplace=True)

othernejm_dup=cnvs[(cnvs.NEJM_Name!='.') & (cnvs.Type=='DUP')]['Sample'].to_list()
df['NEJM_dup_carrier']='.'
df.loc[df.Sample.isin(othernejm_dup), 'NEJM_dup_carrier']='X'
# Add in specific CNV names
nejmdup_dict={}
for s in othernejm_dup:
	names=[i for i in list(cnvs[(cnvs.Sample==s) & (cnvs.Type=='DUP')]['NEJM_Name'].unique()) if i!='.']
	nejmdup_dict[s]=';'.join(names)
df['NEJM_dup_name']=df.Sample.map(nejmdup_dict)
df.fillna('.', inplace=True)

print(df)

# Look for other large rare CNVs
other_del=cnvs[(cnvs.NEJM_Name=='.') & (cnvs.Type=='DEL') & (cnvs.Length>=500000)]['Sample'].to_list()
df['LargeRareCNV_del_carrier']='.'
df.loc[df.Sample.isin(other_del), 'LargeRareCNV_del_carrier']='X'
other_dup=cnvs[(cnvs.NEJM_Name=='.') & (cnvs.Type=='DUP') & (cnvs.Length>=500000)]['Sample'].to_list()
df['LargeRareCNV_dup_carrier']='.'
df.loc[df.Sample.isin(other_dup), 'LargeRareCNV_dup_carrier']='X'

# Save to file
df.to_csv('Intermediate_files/1_cnv_samples.csv', index=False)
