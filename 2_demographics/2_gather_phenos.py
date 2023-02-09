import pandas as pd
import numpy as np

# We want to examine additional phenotypes, so gather those from the raw phenotype files
# Phenotypes:
# BMI (21001)
# Indices of Multiple Depravation (many - https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=76)
# Age at death: Rather than using the age at death from the old phenotype file, I will use dates of birth and dates of death downloaded from the UK Biobank Data Portal (see script 2.1)
field_codes={'eid':'eid',
	'52':'Month of birth',
	'34':'Year of birth',
	'21001':'BMI',
	'26433':'Access to services score (Scotland)',
	'26422':'Access to services score (Wales)',
	'26425':'Community safety score (Wales)',
	'26416':'Crime score (England)',
	'26434':'Crime score (Scotland)',
	'26414':'Education score (England)',
	'26431':'Education score (Scotland)',
	'26421':'Education score (Wales)',
	'26412':'Employment score (England)',
	'26429':'Employment score (Scotland)',
	'26419':'Employment score (Wales)',
	'26413':'Health score (England)',
	'26430':'Health score (Scotland)',
	'26420':'Health score (Wales)',
	'26415':'Housing score (England)',
	'26432':'Housing score (Scotland)',
	'26423':'Housing score (Wales)',
	'26411':'Income score (England)',
	'26428':'Income score (Scotland)',
	'26418':'Income score (Wales)',
	'26410':'Index of Multiple Deprivation (England)',
	'26427':'Index of Multiple Deprivation (Scotland)',
	'26426':'Index of Multiple Deprivation (Wales)',
	'26417':'Living environment score (England)',
	'26424':'Physical environment score (Wales)'
}

fields=list(field_codes.keys())

# Save space by only loading in relevant columns
file=open('/data4/UK_Biobank/download/ukb30075.csv', 'r', encoding='unicode_escape')
for line in file:
	header = line.rstrip().split(',')
	break
file.close()

cols=[]
for f in fields:
	cols += [i.replace('"', '') for i in header if (f+'-')==i[1:len(f)+2]]

phenos = pd.read_csv('/data4/UK_Biobank/download/ukb30075.csv', encoding='unicode_escape', usecols = ['eid']+cols)
print(phenos.shape)

# Restrict to only case/control samples
df=pd.read_csv('../0_get_carriers/UKB_participant_table.csv')
phenos=phenos[phenos.eid.isin(df.Sample.to_list())]
phenos.dropna(axis=1, how='all', inplace=True)
phenos.dropna(axis=0, thresh=2, inplace=True)
phenos.dropna(axis=1, how='all', inplace=True)

print(phenos.shape)

# Consolodate data
for key in field_codes.keys():
	count=len([i for i in cols if key==i[0:len(key)+1]])
	if count>1:
		print(field_codes[key], count)
		print([i for i in cols if key in i])

# BMI has  multiple columns
# For BMI, use most recent data
def bmi(row):
	bmis=[i for i in row.to_list() if i==i]
	if len(bmis)==0:
		return np.nan
	return bmis[-1]
bmi_cols=[i for i in phenos.columns.to_list() if '21001' in i]
phenos['21001-merge']=phenos[bmi_cols].apply(bmi, axis=1)

phenos.drop(bmi_cols, axis=1, inplace=True)

cols=phenos.columns.to_list()
row_dict={}
for key in field_codes:
	if key!='eid':
		col=[i for i in cols if key+'-'==i[0:len(key)+1]]
		row_dict[col[0]]=field_codes[key]
	else:
		row_dict['eid']='eid'
phenos.columns=[row_dict[i] for i in cols]

# Add date of death
dod=pd.read_csv('analysis_files/death_query.txt', sep='\t')
dod=dod[dod.eid.isin(phenos.eid.to_list())]
dod_dict=dict(zip(dod.eid.to_list(), dod.date_of_death.to_list()))
phenos['DOD']=phenos.eid.map(dod_dict)

# Calculate current age and age at death
import datetime as datetime
phenos['DOB']='15/'+phenos['Month of birth'].astype(str)+'/'+phenos['Year of birth'].astype(str)
dobs=pd.to_datetime(phenos['DOB'], dayfirst=True)
dods=pd.to_datetime(phenos['DOD'], yearfirst=True)
phenos['current_age']=dobs.apply(lambda dob: (datetime.datetime(2023, 2, 8)-dob).days)

df2=pd.DataFrame({'DOB':dobs, 'DOD':dods})
phenos['age_at_death']=df2.apply(lambda x: (x.DOD-x.DOB).days, axis=1)

# Add the cause of death
# Some samples have multiple causes listed - report both separated by :
cause_dict={}
for eid in list(dod.eid.unique()):
	causes=':'.join(dod[dod.eid==eid]['cause_icd10'].to_list())
	cause_dict[eid]=causes
phenos['cause_of_death']=phenos.eid.map(cause_dict)

# Annotate case/control status
cc_dict=dict(zip(df.Sample.to_list(), df.Case_Control.to_list()))
phenos['Case_Control']=phenos.eid.map(cc_dict)

# Get counts of each phenotype within each case/control group
controls=['NoCNV_Control', 'LargeRare_Control', 'NEJM_Control']
for c in ['Case']+controls:
	print(phenos[phenos.Case_Control==c].count())

# Due to limited sample size in the cases, I'll drop the indices of depravation in Scotland (n=15) and Wales (n=9)
cols=[i for i in phenos.columns.to_list() if 'Wales' not in i and 'Scotland' not in i]
phenos=phenos[cols]
print(phenos)

# Save to file
phenos.to_csv('analysis_files/2_additional_phenotypes.csv', index=False)
