import pandas as pd

# Combine the tables separated by chromosome into one table
chroms=pd.read_csv('analysis_files/chromosomes.list', header=None, names=['chroms']).chroms.to_list()

outdf=pd.DataFrame({'Sample':[], 'Sex':[], 'Case_Control':[]})
gene_chromosome={}
for c in chroms:
	# Skip chromosome Y
	if c=='chrY':
		continue
	df=pd.read_csv('analysis_files/3_secondhit_tables/'+c+'_secondhit_table.csv')
	gene_chromosome[c] = df.columns.to_list()
	outdf=pd.merge(outdf, df, on=['Sample', 'Sex', 'Case_Control'], how='outer')
	print(c)
	print(outdf)

outdf.fillna(0, inplace=True)

# Save to file
# To make files smaller and more managable, save each chromosome separately
for c in chroms:
	# Skip chromosome Y
	if c=='chrY':
		continue
	outdf[gene_chromosome[c]].to_csv('analysis_files/4_final_secondhit_tables/4_'+c+'_secondhit_table.csv', index=False)
	print(c)

