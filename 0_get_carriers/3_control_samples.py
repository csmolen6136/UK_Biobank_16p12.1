import pandas as pd

# Add case/control labels to the table
df=pd.read_csv('Intermediate_files/2_phenotype_samples.csv')

# Remove missing samples (those with negative sample IDs)
df=df[df.Sample>=0]

# Designate cases (16p12.1 del carriers) and controls
# We have 3 kinds of controls:
# Those with other NEJM deletions
# Those with other large rare deletions
# And those with no large CNVs
df['Case_Control']='.'
df.loc[(df.CNV_data_available=='X') & (df.LargeRareCNV_dup_carrier=='.') & (df.NEJM_dup_carrier=='.'), 'Case_Control']='NoCNV_Control'
df.loc[df.LargeRareCNV_del_carrier=='X', 'Case_Control']='LargeRare_Control'
df.loc[df.NEJM_del_carrier=='X', 'Case_Control']='NEJM_Control'
df.loc[df['16p12.1del_carrier']=='X', 'Case_Control']='Case'

# For second hit analysis, we only care about cases and controls
# Make a list of cases and controls for variant annotation
case_control=df[df.Case_Control!='.'][['Sample', 'Case_Control']]
# Save to file
case_control.to_csv('Intermediate_files/3_case_control.list', index=False, sep='\t')

# Save full table to file
df.to_csv('Intermediate_files/3_control_samples.csv', index=False)
