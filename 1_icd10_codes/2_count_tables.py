import pandas as pd

# Create a table of counts for each node, Sub-block, Block, and Chapter
icd10=pd.read_csv('analysis_files/1_icd_meaning.csv')

df=pd.read_csv('../0_get_carriers/UKB_participant_table.csv')
df=df[df.ICD10=='X']

type_ids={}
types=['Case', 'NoCNV_Control', 'LargeRare_Control', 'NEJM_Control']
for type in types:
	type_ids[type]=df[df.Case_Control==type]['Sample'].to_list()

icd10['Case_Control']=icd10.Sample.map(dict(zip(df.Sample.to_list(), df.Case_Control.to_list())))

icd10.sort_values(by=['Chapter', 'Block', 'Sub-block', 'node'], inplace=True)
print(icd10)

# Get case/control counts
def get_counts(id, category):
	ids_with=icd10[icd10[category]==id]['Sample'].to_list()
	out=[id, category]
	for type in types:
		type_with=list(set(type_ids[type]).intersection(set(ids_with)))
		type_wo=list(set(type_ids[type])-set(ids_with))
		out+=[len(type_with), len(type_wo)]
	return out

out_lst=[]
# Nodes
for node in list(icd10.node.unique()):
	if node=='.':
		continue
	out_lst.append(get_counts(node, 'node'))

# Sub-blocks
for sb in list(icd10['Sub-block'].unique()):
	if sb=='.':
		continue
	out_lst.append(get_counts(sb, 'Sub-block'))
# Blocks
for block in list(icd10.Block.unique()):
	if block=='.':
		continue
	out_lst.append(get_counts(block, 'Block'))
# Chapter
for ch in list(icd10.Chapter.unique()):
	out_lst.append(get_counts(ch, 'Chapter'))

out_df=pd.DataFrame(out_lst, columns=['ID', 'category', 'Case_w', 'Case_wo', 'NoCNV_Control_w', 'NoCNV_Control_wo', 'LargeRare_Control_w', 'LargeRare_Control_wo', 'NEJM_Control_w', 'NEJM_Control_wo'])
# Save to file
out_df.to_csv('analysis_files/2_count_table.csv', index=False)
