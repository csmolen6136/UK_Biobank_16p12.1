import pandas as pd
import numpy as np

# Add in CNV second hits (as number of genes del or dup)
df=pd.read_csv('Intermediate_files/4_snv_samples.csv')
print(df)

cnvs=pd.read_csv('/data5/UK_Biobank/PennCNV_calls/bed_files/19_gene_filter.bed', sep='\t')
cnvs['Num_genes']=cnvs.gene_ids.str.count(' ')+1
print(cnvs)

def count_genes(row, type='DEL'):
	# Skip samples with no CNV data
	if row['CNV_data_available']!='X':
		return np.nan
	# Skip samples that will not be used
	cc=row['Case_Control']
	if cc=='.':
		return np.nan
	sample=row['Sample']
	subdf=cnvs[(cnvs.Sample==sample) & (cnvs.Type==type)]
	if type=='DEL':
		subdf=subdf[(subdf.NEJM_Name!='16p12.1')]
		if cc=='NEJM_Control':
			# If they have multiple NEJM CNVs, just pick the first as the first hit
			nejm_cnvs=row['NEJM_del_name'].split(';')
			subdf=subdf[subdf.NEJM_Name!=nejm_cnvs[0]]
		if cc=="LargeRare_Control":
			# If they have multiple large, rare CNVs, pick the largest as the first hit
			maxlen=max(subdf[subdf.Length>=500000]['Length'].to_list())
			subdf=subdf[subdf.Length!=maxlen]
	return sum(subdf.Num_genes.to_list())

df['Genes_del']=df.apply(count_genes, axis=1)
df['Genes_dup']=df.apply(count_genes, type='DUP', axis=1)

print(df)

# Save to file
df.to_csv('Intermediate_files/5_cnv_secondhit.csv', index=False)
