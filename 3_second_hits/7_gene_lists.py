import pandas as pd

# Create lists of second-hit genes for each case-control group and save to file
df=pd.read_csv('analysis_files/2_secondhits.csv')

controls=['NoCNV_Control', 'LargeRare_Control', 'NEJM_Control']

for c in ['Case']+controls:
	subdf=df[df.Case_Control==c]

	for sex in ['M', 'F']:
		genes=list(set(df[(df.Case_Control==c) & (df.Sex==sex)]['gene_id'].to_list()))
		genes2=[i for i in genes if ';' not in i]
		genes3=[i for i in genes if ';' in i]
		for g in genes3:
			genes2+=';'.split(g)
		genes2.sort()

		outfile=open('analysis_files/7_gene_lists/'+c+'_'+sex+'.list', 'w')
		outfile.write('\n'.join(['gene_id']+genes2))
		outfile.close()

