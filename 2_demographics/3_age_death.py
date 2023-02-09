import pandas as pd
import numpy as np

from lifelines import KaplanMeierFitter

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Look for differences in age at death between cases and controls
df=pd.read_csv('analysis_files/2_additional_phenotypes.csv')
df=df[['eid', 'current_age', 'age_at_death', 'Case_Control']]

df['duration']=df.age_at_death
df.loc[df.age_at_death.isnull(), 'duration']=df.current_age
# Convert duration from days to years
df.duration=df.duration/365.25

df['status']=1
df.loc[df.age_at_death.isnull(), 'status']=0

df.sort_values(by=['Case_Control', 'duration'], inplace=True)

print(df)

# Estimate Kaplan-Meier survivability curve for each case/control group
controls=['NoCNV_Control', 'LargeRare_Control', 'NEJM_Control']

pdf=PdfPages('figures/3_KM_survival_curves.pdf')

for cont in controls:
	for c in ['Case', cont]:
		subdf=df[df.Case_Control==c]

		# Survival Curve
		kmf = KaplanMeierFitter()
		kmf.fit(subdf.duration, subdf.status)

		survival=kmf.survival_function_
		ci=kmf.confidence_interval_survival_function_

		out=survival.join(ci)

		# Save to file
		out.to_csv('result_tables/3_KM_survival_'+c+'.csv')

		# Plot curve and CI
		if c=='Case':
			color='C0'
		else:
			color='C1'
		ages=out.index
		lo, hi = np.transpose(ci.values)

		plt.fill_between(ages, lo, hi, color=color, alpha=0.2)
		plt.plot(ages, out.KM_estimate, color=color, label=c)

	plt.ylabel('Survival')
	plt.xlabel('Age (years)')
	plt.legend()

	pdf.savefig()
	plt.close()

pdf.close()
