import pandas as pd

# Tables of ICD codes were generated as part of creating the participants table
# Take tables of codes and annotate them with meanings
icd10=pd.read_csv('../0_get_carriers/Phenotype_files/2_icd10_codes.csv')

# This file was made before the list of participants was finalized
# To save time, limit to only individuals relevant for this analysis
df=pd.read_csv('../0_get_carriers/UKB_participant_table.csv')
icd10=icd10[icd10.Sample.isin(df.Sample.to_list())]

# Meaning was downloaded from: https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=19
coding=pd.read_csv('analysis_files/coding19.tsv', sep='\t')

# Add specific meaning
icd10['meaning']=icd10.ICD10_code.map(dict(zip(coding.coding.to_list(), coding.meaning.to_list())))

# Add node ID
icd10['node']=icd10.ICD10_code.map(dict(zip(coding.coding.to_list(), coding.node_id.to_list())))

# Add categories to the meanings - 3 labels
# Chapter (i.e. "Chapter I Certain infectious and parasitic diseases")
# Block (i.e. "A00-A09 Intestinal infectious diseases")
# Sub-block (i.e. "A00 Cholera")
coding['label']='.'
coding.loc[coding.coding.str.contains('Chapter'), 'label']='Chapter'
coding.loc[coding.coding.str.contains('Block'), 'label']='Block'
coding.loc[~(coding.meaning.str.contains('\\.')) & (~coding.coding.str.contains('Block')) & (~coding.coding.str.contains('Chapter')), 'label']='Sub-block'
print(coding.label.value_counts())

# Assign each code to a Chapter, Block, and Sub-block
def get_label(row, category='Chapter'):
	node=row['node_id']
	parent=row['parent_id']
	label=row['label']
	if label==category:
		return node
	if label=='Chapter' and category!='Chapter':
		return '.'
	if label=='Block' and category=='Sub-block':
		return '.' 
	else:
		row_parent=coding[coding.node_id==parent][['node_id', 'parent_id', 'label']].squeeze()
		return get_label(row_parent, category=category)
coding['Chapter']=coding[['node_id', 'parent_id', 'label']].apply(get_label, category='Chapter', axis=1)
print(coding)
coding['Block']=coding[['node_id', 'parent_id', 'label']].apply(get_label, category='Block', axis=1)
print(coding)
coding['Sub-block']=coding[['node_id', 'parent_id', 'label']].apply(get_label, category='Sub-block', axis=1)
print(coding)

# Save this structure to file
coding.to_csv('analysis_files/coding19_updated.csv', index=False)

# Add Chapter, Block, and Sub-block to ICD10 codes
icd10['Sub-block']=icd10.ICD10_code.map(dict(zip(coding.coding.to_list(), coding['Sub-block'].to_list())))
icd10['Block']=icd10.ICD10_code.map(dict(zip(coding.coding.to_list(), coding['Block'].to_list())))
icd10['Chapter']=icd10.ICD10_code.map(dict(zip(coding.coding.to_list(), coding['Chapter'].to_list())))
print(icd10)

# Save to file
icd10.to_csv('analysis_files/1_icd_meaning.csv', index=False)
