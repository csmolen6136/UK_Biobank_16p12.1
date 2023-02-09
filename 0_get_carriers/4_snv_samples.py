import pandas as pd

# Create a table of samples in the UK Biobank
# Table contains ID, sex, 16p12.1 carrier status, other NEJM CNV carrier status, other large (>500 kb) rare genic CNV carrier status, ICD10 code availability, 2nd hits availability
# This script will only add in SNV data

df=pd.read_csv('Intermediate_files/3_control_samples.csv')

burden=pd.read_csv('/data5/UK_Biobank/annotations/annovar/2023_01_05/tables/19_burden_table.csv', header=None, index_col=0, names=['burden'])

df['SNV_data_available']='.'
df.loc[df.Sample.isin(burden.index.to_list()), 'SNV_data_available']='X'

df['SNV_burden']=df.Sample.map(burden.burden.to_dict())

print(df[df.SNV_data_available=='X'])

# Save to file
df.to_csv('Intermediate_files/4_snv_samples.csv', index=False)
