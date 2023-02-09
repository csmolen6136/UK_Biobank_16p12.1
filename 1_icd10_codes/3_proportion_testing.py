import pandas as pd
from statsmodels.stats.proportion import proportions_ztest

# Compare the proportions of cases and controls with and without each ICD10 ID
df=pd.read_csv('analysis_files/2_count_table.csv')

controls=['NoCNV_Control', 'LargeRare_Control', 'NEJM_Control']

out=[]
for idx, row in df.iterrows():
	case_w=row['Case_w']
	case_wo=row['Case_wo']

	for cont in controls:
		cont_w=row[cont+'_w']
		cont_wo=row[cont+'_wo']

		statistic, p = proportions_ztest([case_w, cont_w], [(case_w+case_wo), (cont_w+cont_wo)], alternative='two-sided', value=0)
		diff=(case_w/(case_w+case_wo))-(cont_w/(cont_w+cont_wo))

		out.append([row['ID'], row['category'], cont, case_w, case_wo, cont_w, cont_wo, diff, statistic, p])

output=pd.DataFrame(out, columns=['ID', 'category', 'Control', 'Case_w', 'Case_wo', 'Control_w', 'Control_wo', 'difference', 'statistic', 'p'])

# Multiple test Correct by control and category
factor_dict={}
for cont in controls:
	for type in ['node', 'Sub-Block', 'Block', 'Chapter']:
		factor_dict[cont+'.'+type]=output[(output.category==type) & (output.Control==cont) & ~(output.statistic.isnull())].shape[0]
output['type_cont']=output.Control+'.'+output.category
output['factor']=output.type_cont.map(factor_dict)
output['bonferroni_p']=output.p*output.factor
output.loc[output.bonferroni_p>1, 'bonferroni_p']=1

output=output[['ID', 'category', 'Control', 'Case_w', 'Case_wo', 'Control_w', 'Control_wo', 'difference', 'statistic', 'p', 'bonferroni_p']]

# Add in whether term is enriched or depleted
output['direction']='enriched'
output.loc[output.difference<=0, 'direction']='depleted'

# Add in meaning
coding=pd.read_csv('analysis_files/coding19_updated.csv')
meaning_dict=dict(zip(coding.node_id.to_list(), coding.meaning.to_list()))
output['ID']=output.ID.astype(int)
output['meaning']=output.ID.map(meaning_dict)

print(output[output.bonferroni_p<=0.05])
print(output[output.bonferroni_p<=0.05]['direction'].value_counts())
print(output[output.bonferroni_p<=0.05]['Control'].value_counts())
print(output[output.bonferroni_p<=0.05]['category'].value_counts())

# Save to file
output.to_csv('result_tables/3_proprtion_z_results.csv', index=False)
