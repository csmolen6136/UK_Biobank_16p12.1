import pandas as pd

# Finalize participant table

df=pd.read_csv('Intermediate_files/5_cnv_secondhit.csv')

# Restrict only to cases and controls
df=df[df.Case_Control!='.']

# To be useful in the analysis, samples need to have phenotype or SNV and CNV data
df=df[((df.SNV_data_available=='X') & (df.CNV_data_available)) | (df.ICD10=='X')]

# Reformat some columns
df.Genes_del=df.Genes_del.astype(int)
df.Genes_dup=df.Genes_dup.astype(int)

print(df.Case_Control.value_counts())

# Limit to relevant columns
df=df[['Sample', 'Sex', '16p12.1del_carrier', 'Case_Control', 'CNV_data_available', 'SNV_data_available', 'ICD10', 'NEJM_del_name', 'LargeRareCNV_del_carrier', 'SNV_burden', 'Genes_del', 'Genes_dup']]

# Save to file
df.to_csv('UKB_participant_table.csv', index=False)
