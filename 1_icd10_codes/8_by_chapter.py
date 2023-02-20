import pandas as pd

# Re-calculate FDR values while only considering nodes, blocks, and sub-blocks within significant chapters
# Consider only one chapter at a time
df=pd.read_csv('result_tables/6_parsed_model.csv')

# Add in chapter information
coding=pd.read_csv('analysis_files/coding19_updated.csv')
df['chapter']=df.id.map(dict(zip(coding.node_id.to_list(), coding.Chapter.to_list())))

# Restrict dataframe to only IDs within significant chapters
controls=['NoCNV_Control', 'LargeRare_Control', 'NEJM_Control']
control_dict={}
for c in controls:
	control_dict[c]=df[(df.control==c) & (df.type=='Chapter') & (df.cc_bonferroni_fdr<=0.05)]['id'].to_list()
df['keep']=False
for c in controls:
	df.loc[(df.control==c) & (df.chapter.isin(control_dict[c])), 'keep']=True
df=df[df.keep]

# Re-calculate FDR considering one chapter-control-type combination at a time
chaps=list(df.chapter.unique())
types=['node', 'Sub-block', 'Block', 'Chapter']
df['factor']=0
for c in controls:
	for ch in chaps:
		for t in types:
			count=df[(df.control==c) & (df.chapter==ch) & (df.type==t)].shape[0]
			df.loc[(df.control==c) & (df.chapter==ch) & (df.type==t), 'factor']=count*2

print(df[df.cc_bonferroni_fdr<=0.05][['type', 'id', 'control', 'cc_bonferroni_fdr', 'cc_posneg', 'sex_bonferroni_fdr', 'sex_posneg', 'meaning']])

df['cc_bonferroni_fdr']=df.cc_p*df.factor
df.loc[df.cc_bonferroni_fdr>1, 'cc_bonferroni_fdr']=1
df['sex_bonferroni_fdr']=df.sex_p*df.factor
df.loc[df.sex_bonferroni_fdr>1, 'sex_bonferroni_fdr']=1

df.sort_values(by=['control', 'chapter', 'type'], inplace=True)

print(df[df.cc_bonferroni_fdr<=0.05][['type', 'id', 'control', 'cc_bonferroni_fdr', 'cc_posneg', 'sex_bonferroni_fdr', 'sex_posneg', 'meaning']])

# Save to file
df.to_csv('result_tables/8_by_chapter.csv', index=False)

# Plot significant results
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Plot significant results within each control-chapter
df=df[df.cc_bonferroni_fdr<=0.05]
df['type_cat']=pd.Categorical(df['type'], categories=types, ordered=True)
df.sort_values(by=['type_cat', 'meaning'], inplace=True)

pdf=PdfPages('figures/8_by_chapter.pdf')
for c in controls:
	for ch in chaps:
		subdf=df[(df.control==c) & (df.chapter==ch)]
		if subdf.shape[0]==0:
			continue
		# Make a forest plot showing all significant results
		ys=list(range(0, subdf.shape[0]))
		# Plot Case-Control estimates and error bars
		plt.scatter(subdf.cc_est.to_list(), [i+0.1 for i in ys], color='#2c3e50', label='Case-Control')
		for i in ys:
			est=subdf.cc_est.to_list()[i]
			err=subdf.cc_err.to_list()[i]
			plt.plot([est-err, est+err], [i+0.1, i+0.1], color='#2c3e50')
		# Plot sex estimates and error bars
		# Some codes have very high error values for sex - maintain view of current plot
		lo, hi = plt.xlim()
		plt.scatter(subdf.sex_est.to_list(), [i-0.1 for i in ys], color='#f1c40f', label='Sex')
		for i in ys:
			est=subdf.sex_est.to_list()[i]
			err=subdf.sex_err.to_list()[i]
			plt.plot([est-err, est+err], [i-0.1, i-0.1], color='#f1c40f')
		if max(subdf.sex_err.to_list())>3:
			plt.xlim(lo-2, hi+2)
		# Add axis label
		meanings=subdf.meaning.to_list()
		ylabels=[]
		for i in meanings:
			if len(i)>60:
				ylabels.append(i[0:60]+' . . .')
			else:
				ylabels.append(i)
		plt.yticks(ys, ylabels)
		# Add title
		plt.title(c+' '+str(ch))
		# Add a line at 0
		plt.plot([0, 0], [-0.5, max(ys)+0.5], color='k', ls=':', zorder=0)
		plt.ylim(-0.25, max(ys)+0.25)

		plt.tight_layout()
		pdf.savefig()
		plt.close()

pdf.close()
