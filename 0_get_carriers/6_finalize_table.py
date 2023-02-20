import pandas as pd

# Finalize participant table

df=pd.read_csv('Intermediate_files/5_cnv_secondhit.csv')

# Restrict only to cases and controls
df=df[df.Case_Control!='.']

# To be useful in the analysis, samples need to have phenotype or SNV and CNV data
df=df[((df.SNV_data_available=='X') & (df.CNV_data_available)) | (df.ICD10=='X')]

# Add ethnicity data
phenos = pd.read_csv('/data4/UK_Biobank/download/ukb30075.csv', encoding='unicode_escape', usecols = ['eid', '21000-0.0', '21000-1.0', '21000-2.0'])
phenos=phenos[phenos.eid.isin(df.Sample.to_list())]
def merge_codes(visits):
	codes=list(set([i for i in visits if i==i and i>0]))
	if len(codes)==1:
		return codes[0]
	# Codes can be hierarchical, if codes are not conflicting in the same hierarchy, return most specific code
	first_nums = list(set([str(i)[0] for i in codes]))
	if len(first_nums)==1 and float(first_nums[0]) in codes and len(codes)==2:
		return max(codes)
	# If more than one code or conflicting hierarchical codes are reported, return unknown
	return -1
phenos['code_merge']=phenos[['21000-0.0', '21000-1.0', '21000-2.0']].apply(merge_codes, axis=1)
# Map to ethnicity using coding 1001
coding=pd.read_csv('Intermediate_files/coding1001.tsv', sep='\t')
phenos['ethnicity']=phenos.code_merge.map(dict(zip(coding.coding.to_list(), coding.meaning.to_list())))

df['ethnicity']=df.Sample.map(dict(zip(phenos.eid.to_list(), phenos.ethnicity.to_list())))
print(df.ethnicity.value_counts())

# Reformat some columns
df.Genes_del=df.Genes_del.astype(int)
df.Genes_dup=df.Genes_dup.astype(int)

print(df.Case_Control.value_counts())

# Limit to relevant columns
df=df[['Sample', 'Sex', 'ethnicity', '16p12.1del_carrier', 'Case_Control',
	'CNV_data_available', 'SNV_data_available', 'ICD10',
	'NEJM_del_name', 'LargeRareCNV_del_carrier',
	'SNV_burden', 'Genes_del', 'Genes_dup']]

# Save to file
df.to_csv('UKB_participant_table.csv', index=False)
