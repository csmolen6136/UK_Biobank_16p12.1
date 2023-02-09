import pandas as pd

# Gather specific second hits for further analysis
df=pd.read_csv('../0_get_carriers/UKB_participant_table.csv')

# Restrict to only individuals with both CNV and SNV data
df=df[df.SNV_data_available=='X']
cc_dict=dict(zip(df.Sample.to_list(), df.Case_Control.to_list()))

# Get CNV second hits
cnvs=pd.read_csv('/data5/UK_Biobank/PennCNV_calls/bed_files/19_gene_filter.bed', sep='\t')
cnvs=cnvs[cnvs.Sample.isin(df.Sample.to_list())]
# Remove first hits
cnvs=cnvs[~((cnvs.NEJM_Name=='16p12.1') & (cnvs.Type=='DEL'))]
cnvs['Case_Control']=cnvs.Sample.map(cc_dict)

def anno_firsthit(row):
	type=row['Type']
	if type=='DUP':
		return '.'

	cc=row['Case_Control']
	if cc=='Case' or cc=='NoCNV_Control':
		return '.'

	sample=row['Sample']
	if cc=='NEJM_Control':
		nejm=df[df.Sample==sample]['NEJM_del_name'].to_list()[0].split(';')
		# If they have multiple NEJM CNVs, the first is the first hit
		if row['NEJM_Name']==nejm[0]:
			return 'X'
		return '.'

	if cc=='LargeRare_Control':
		# If they have multiple large, rare CNVs, pick the largest as the first hit
		maxlen=max(cnvs[(cnvs.Sample==sample) & (cnvs.Length>=500000) & (cnvs.Type=='DEL')]['Length'].to_list())
		if row['Length']==maxlen:
			return 'X'
		return '.'

cnvs['firsthit']=cnvs.apply(anno_firsthit, axis=1)
print(cnvs.firsthit.value_counts())

cnvs=cnvs[cnvs.firsthit=='.']

# Reformat CNV data on a gene by gene level
cnvs['Ref']='CN2'
cnvs['Alt']=cnvs.Type
cnvs['Mut_type']='CNV_'+cnvs.Type
cnvs['frequency']=cnvs.UKB_freq

var_lst=[]
for idx, row in cnvs.iterrows():
	gene_ids=row['gene_ids'].split(' ')
	gene_names=row['gene_symbols'].split(' ')

	for g in range(len(gene_ids)):
		var_lst.append([row['Sample'], row['Case_Control'], row['Chr'], row['Start'], row['End'], row['Ref'], row['Alt'], row['variant_id'],
				row['Mut_type'], gene_ids[g], gene_names[g], row['frequency']])

var_df=pd.DataFrame(var_lst, columns=['Sample', 'Case_Control', 'Chr', 'Start', 'End', 'Ref', 'Alt', 'variant_id', 'Mut_type', 'gene_id', 'gene_symbol', 'frequency'])


# Add in SNV data
snv=pd.read_csv('/data5/UK_Biobank/annotations/annovar/2023_01_05/tables/18_gene_ids.csv')

def remove_dot(id):
	if ';' not in id:
		return id.split('.')[0]

	ids=id.split(';')
	ids_out=[]
	for i in ids:
		ids_out.append(i.split('.')[0])

	return ';'.join(ids_out)

snv['gene_id']=snv.Gene_id.apply(remove_dot)
snv['gene_symbol']=snv.Gene_symbol
snv['Start']=snv.Pos
snv['End']=snv.Pos
snv['Case_Control']=snv.Sample.map(cc_dict)

var_df=pd.concat([var_df, snv[['Sample', 'Case_Control', 'Chr', 'Start', 'End', 'Ref', 'Alt', 'variant_id', 'Mut_type', 'gene_id', 'gene_symbol', 'frequency']]])

var_df.sort_values(by=['Case_Control', 'Sample', 'Chr', 'Start'], inplace=True)

# Add sex
sex_dict=dict(zip(df.Sample.to_list(), df.Sex.to_list()))
var_df['Sex']=var_df.Sample.map(sex_dict)

var_df=var_df[['Sample', 'Case_Control', 'Sex', 'Chr', 'Start', 'End', 'Ref', 'Alt', 'variant_id', 'Mut_type', 'gene_id', 'gene_symbol', 'frequency']]

# Save to file
var_df.to_csv('analysis_files/2_secondhits.csv', index=False)
