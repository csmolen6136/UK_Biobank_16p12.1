import pandas as pd
import numpy as np

# Gather the cognitive phenotypes used in the Kendall paper (https://pubmed.ncbi.nlm.nih.gov/27773354/)

# Kendall
# 20132 - Errors in pairs matching test
# 20023 - Mean time to correctly identify matches (reaction time)
# 20016 - Fluid intelligence score
# 20240 - Maximum digits remembered correctly (numeric memory)
# 20159 - Number of symbol-digit matches made correctly (symbol digit substitution - complex processing speed)
# 20156 - Time to complete test (Trail making A - visual attention)
# 20157 - Time to complete test (Trail making B - visual attention)

cognitive={'eid':'Sample',
		'20132':'pairs_matching_errors',
		'20023':'reaction_time_mean',
		'20016':'fluid_intelligence',
		'20240':'numeric_memory_maxdigits',
		'20159':'sd_sub_correctmatches',
		'20156':'traila_time',
		'20157':'trailb_time'}

# Save space by only loading in relevant columns
file=open('/data4/UK_Biobank/download/ukb30075.csv', 'r', encoding='unicode_escape')
for line in file:
	header = line.rstrip().split(',')
	break
file.close()

cols=[]
for f in list(cognitive.keys()):
	cols += [i.replace('"', '') for i in header if (f+'-')==i[1:len(f)+2]]

phenos = pd.read_csv('/data4/UK_Biobank/download/ukb30075.csv', encoding='unicode_escape', usecols = ['eid']+cols)

# Restrict to only case/control samples
df=pd.read_csv('../0_get_carriers/UKB_participant_table.csv')
phenos=phenos[phenos.eid.isin(df.Sample.to_list())]
phenos.dropna(axis=1, how='all', inplace=True)
phenos.dropna(axis=0, thresh=2, inplace=True)
phenos.dropna(axis=1, how='all', inplace=True)

print(phenos)
for key in cognitive.keys():
	count=len([i for i in cols if key in i])
	if count>1:
		print(cognitive[key], count)
		print([i for i in cols if key in i])

# Pairs matching, reaction time, and fluid intelligence all have multiple values
# For fluid intelligence I will only take the first instance, as it has a different distribution than other instances (https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20016)
# For reaction time, I will only take the first instance for the same reason (https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20023)
# For pairs matching, I will take the first score if there are multiple scores
drop_cols=[i for i in cols if ('20016' in i and i!='20016-0.0') or ('20023' in i and i!='20023-0.0')]
phenos.drop(drop_cols, axis=1, inplace=True)

def pairs_matching_merge(merge_cols):
	scores=[i for i in merge_cols.to_list() if i==i]
	if len(scores)>0:
		return scores[0]
	return np.nan

pm_cols=[i for i in cols if '20132' in i]
phenos['20132_merge']=phenos[pm_cols].apply(pairs_matching_merge, axis=1)
phenos.drop(pm_cols, axis=1, inplace=True)

print(phenos)
cols=phenos.columns.to_list()
row_dict={}
for key in cognitive.keys():
	col=[i for i in cols if key in i]
	row_dict[col[0]]=cognitive[key]
phenos.columns=[row_dict[i] for i in cols]

# Add case/control status
cc_dict=dict(zip(df.Sample.to_list(), df.Case_Control.to_list()))
phenos['Case_Control']=phenos.Sample.map(cc_dict)

# Get counts of each phenotype within each case/control group
controls=['NoCNV_Control', 'LargeRare_Control', 'NEJM_Control']
for c in ['Case']+controls:
	print(phenos[phenos.Case_Control==c].count())

# Save to file
phenos.to_csv('analysis_files/1_kendall_cognitive.csv', index=False)
