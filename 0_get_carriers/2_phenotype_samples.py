import pandas as pd

# Create a table of samples in the UK Biobank
# Table contains ID, sex, 16p12.1 carrier status, other NEJM CNV carrier status, other large (>500 kb) rare genic CNV carrier status, ICD10 code availability, 2nd hits availability
# This script will add in phenotype data

df=pd.read_csv('Intermediate_files/1_cnv_samples.csv')

# To save space, only load in relevant columns of phenotype file
file=open('/data4/UK_Biobank/download/ukb30075.csv', 'r', encoding='unicode_escape')
for line in file:
	header = line.rstrip().split(',')
	break
file.close()

# There are 4 ICD Data Fields
# 41202 - Main ICD10
# 41203 - Main ICD9
# 41204 - Secondary ICD10
# 41205 - Secondary ICD9
icd=['41202', '41203', '41204', '41205']
cols=[]
for ic in icd:
	cols += [i.replace('"', '') for i in header if ic+'-' in i]

phenos = pd.read_csv('/data4/UK_Biobank/download/ukb30075.csv', encoding='unicode_escape', usecols = cols+['eid', '31-0.0'])

# Annotate with presence/absence of ICD information
icd10_cols=[i for i in cols if '41202' in i or '41204' in i]
icd9_cols=[i for i in cols if '41203' in i or '41205' in i]

phenos['ICD10']='X'
phenos.loc[phenos[icd10_cols].isnull().apply(lambda x: all(x), axis=1), 'ICD10']='.'

phenos['ICD9']='X'
phenos.loc[phenos[icd9_cols].isnull().apply(lambda x: all(x), axis=1), 'ICD9']='.'
print(phenos.ICD9.value_counts())
print(phenos.ICD10.value_counts())

missing=list(set(df.Sample.to_list())-set(phenos.eid.to_list()))
print(len(missing))
print(missing)

# Only samples missing are those tht have been removed from the UK Biobank (negative values) -- good!
phenos['Sample']=phenos['eid']
phenos['Sex']=phenos['31-0.0'].map({0:'F', 1:'M'})

df=pd.merge(df, phenos[['Sample', 'Sex', 'ICD10', 'ICD9']], on=['Sample'], how='outer')
df.fillna('.', inplace=True)
print(df)

# Save table to file
df.sort_values(by=['Sample'], inplace=True)
df.to_csv('Intermediate_files/2_phenotype_samples.csv', index=False)

# Also save ICD codes to a file
# Reformat the codes into a flat file
icd10=phenos[phenos.ICD10=='X'][['Sample']+icd10_cols]
df_lst=[]
for col in icd10_cols:
	subdf=icd10[['Sample', col]]
	subdf.columns=['Sample', 'ICD10_code']
	subdf=subdf[~subdf.ICD10_code.isnull()]
	if '41202' in col:
		subdf['Type']='Main'
	else:
		subdf['Type']='Secondary'
	if subdf.shape[0]>0:
		df_lst.append(subdf.copy())
icd10_df=pd.concat(df_lst)
print(icd10_df)
# Save to file
icd10_df.to_csv('Phenotype_files/2_icd10_codes.csv', index=False)

# Do the same for ICD9
icd9=phenos[phenos.ICD9=='X'][['Sample']+icd9_cols]
df_lst=[]
for col in icd9_cols:
	subdf=icd9[['Sample', col]]
	subdf.columns=['Sample', 'ICD9_code']
	subdf=subdf[~subdf.ICD9_code.isnull()]
	if '41203' in col:
		subdf['Type']='Main'
	else:
		subdf['Type']='Secondary'
	if subdf.shape[0]>0:
		df_lst.append(subdf.copy())
icd9_df=pd.concat(df_lst)
print(icd9_df)
# Save to file
icd9_df.to_csv('Phenotype_files/2_icd9_codes.csv', index=False)
