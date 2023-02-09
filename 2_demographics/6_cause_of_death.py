import pandas as pd
import scipy.stats as stats

# Look for differences in the cause of death between cases and controls
df=pd.read_csv('analysis_files/2_additional_phenotypes.csv')
df=df[['eid', 'cause_of_death', 'Case_Control']]
df=df[~df.cause_of_death.isnull()]

# Map ICD10 meanings
coding=pd.read_csv('../1_icd10_codes/analysis_files/coding19_updated.csv')
meaning_dict=dict(zip(coding.coding.to_list(), coding.meaning.to_list()))
chapter_dict=dict(zip(coding.coding.to_list(), coding.Chapter.to_list()))
chapter_meaning=dict(zip(coding.node_id.to_list(), coding.meaning.to_list()))

df['meaning']=df.cause_of_death.map(meaning_dict)
df['Chapter']=df.cause_of_death.map(chapter_dict)
df['Chapter_meaning']=df.Chapter.map(chapter_meaning)

# Some samples have multiple ICD10 codes
# Update them with chapters as well
def update_multiple(row, return_what='chapter_meaning'):
	causes=list(set(row['cause_of_death'].split(':')))
	meanings=':'.join([meaning_dict[i] for i in causes])
	if return_what=='meaning':
		return meanings
	chapter_lst=list(set([chapter_dict[i] for i in causes]))
	chapters=':'.join([str(i) for i in chapter_lst])
	if return_what=='chapter':
		return chapters
	chapter_meanings=':'.join(list(set([chapter_meaning[i] for i in chapter_lst])))
	return chapter_meanings

df.loc[df.cause_of_death.str.contains(':'), 'meaning']=df.loc[df.cause_of_death.str.contains(':'), ].apply(update_multiple, return_what='meaning', axis=1)
df.loc[df.cause_of_death.str.contains(':'), 'Chapter']=df.loc[df.cause_of_death.str.contains(':'), ].apply(update_multiple, return_what='chapter', axis=1)
df.loc[df.cause_of_death.str.contains(':'), 'Chapter_meaning']=df.loc[df.cause_of_death.str.contains(':'), ].apply(update_multiple, axis=1)

chapters=[i for i in list(df.Chapter_meaning.unique()) if ':' not in i]
chapters.sort()

controls=['NoCNV_Control', 'LargeRare_Control', 'NEJM_Control']
df2_dict={}
for c in ['Case']+controls:
	subdf=df[df.Case_Control==c]
	out=[]
	for ch in chapters:
		out.append(subdf[subdf.Chapter_meaning.str.contains(ch)].shape[0])
	df2_dict[c]=out

# Perform chi-sq to check for differences in Chapters between groups
df2=pd.DataFrame(df2_dict, index=chapters)
# Save to file
df2.to_csv('result_tables/6_causes_of_death.csv')

stat_lst=[]
for c in controls:
	df2[c]=df2['Case'].sum()*df2[c]/df2[c].sum()
	testdf=df2[['Case', c]]
	# Remove rows where case and control are 0
	testdf=testdf[~((testdf.Case==0) & (testdf[c]==0))]
	chisq, p = stats.chisquare(testdf['Case'], testdf[c])
	stat_lst.append([c, 'chi-squared', chisq, p])

stat_df=pd.DataFrame(stat_lst, columns=['control', 'test', 'statistic', 'p'])
stat_df.to_csv('result_tables/6_causes_chisq.csv', index=False)
