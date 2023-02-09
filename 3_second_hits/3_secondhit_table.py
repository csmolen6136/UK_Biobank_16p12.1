import pandas as pd
import sys


# Convert the secondhit table into a flattened table
chr=sys.argv[1]

# Each gene is a column and each sample is a row
df=pd.read_csv('analysis_files/2_secondhits.csv')

# Limit to one chromosome to save time
df=df[df.Chr==chr]

samples=list(df.Sample.unique())
genes=list(df[~df.gene_id.str.contains(';')]['gene_id'].unique())

multi_genes=list(df[df.gene_id.str.contains(';')]['gene_id'].unique())
for g in multi_genes:
	genes+=g.split(';')
genes=list(set(genes))
print(len(genes))

outdf=pd.DataFrame(0, index=samples, columns=genes)
for g in genes:
	gene_samps=df[df.gene_id.str.contains(g)]['Sample'].to_list()
	outdf.loc[gene_samps, g]=1

outdf['Sample']=outdf.index

# Also add in sex and case/control status
ukb=pd.read_csv('../0_get_carriers/UKB_participant_table.csv')
ukb.index=ukb.Sample
ukb.loc[ukb.Sex=='M', 'Sex']=0
ukb.loc[ukb.Sex=='F', 'Sex']=1
ukb.Sex=ukb.Sex.astype(int)
sex_dict=ukb.Sex.to_dict()
cc_dict=ukb.Case_Control.to_dict()

outdf['Sex']=outdf.Sample.map(sex_dict)
outdf['Case_Control']=outdf.Sample.map(cc_dict)

# Save to file
outdf.to_csv('analysis_files/3_secondhit_tables/'+chr+'_secondhit_table.csv', index=False)
