import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as ss

plt.rc('font', size=12)

p = argparse.ArgumentParser()
p.add_argument("--tidy_spectra", required=True,
                    help="""file containing mutation spectra for \
                            each sample in tidy, long-form format""")
p.add_argument("--out", required=True,
                    help="""name of output file""")
args = p.parse_args()

# read in tidy data
tidy_spectra = pd.read_csv(args.tidy_spectra)

# only consider mutation rates (not fractions)
tidy_rates = tidy_spectra.query('estimate_type == "rate"')

cols = ['bxd_strain_conv', 'haplotype_at_qtl', 'estimate']

# subset dataframe into three separate dataframes
# (one for just C>A mutations, one for all but C>A mutations,
# and one for the overall mutation rate including all types)
df_no_ca = tidy_rates.query('base_mut != "C>A"')
df_ca = tidy_rates.query('base_mut == "C>A"')
df_ca = df_ca[cols]
df_ca['rate_type'] = r'C$\to$A'

df_no_ca = df_no_ca.groupby(['bxd_strain_conv', 'haplotype_at_qtl']).sum().reset_index()
df_no_ca = df_no_ca[cols]
df_no_ca['rate_type'] = r'non-C$\to$A'

df_total = tidy_rates.groupby(['bxd_strain_conv', 'haplotype_at_qtl']).sum().reset_index()
df_total = df_total[cols]
df_total['rate_type'] = r'combined'

mean_hap_0 = np.mean(df_ca.query('haplotype_at_qtl == "B"')['estimate'])
mean_hap_1 = np.mean(df_ca.query('haplotype_at_qtl == "D"')['estimate'])
print (ss.ttest_ind(df_ca.query('haplotype_at_qtl == "B"')['estimate'], 
                    df_ca.query('haplotype_at_qtl == "D"')['estimate'],
                    equal_var=False))
print (mean_hap_0, mean_hap_1, mean_hap_1 / mean_hap_0)



# combine the three dataframes, which now each have a unique
# "rate_type" value
combined = pd.concat([df_total, df_no_ca, df_ca])

combined['Haplotype at QTL'] = combined['haplotype_at_qtl']

sns.set_style('ticks')

f, ax = plt.subplots()

sns.boxplot(x="rate_type", y="estimate", hue="Haplotype at QTL", 
                data=combined, ax=ax, color='w', fliersize=0)
sns.stripplot(x="rate_type", y="estimate", hue="Haplotype at QTL", 
                data=combined, ax=ax, palette='colorblind', dodge=True, 
                edgecolor='k', lw=1.5)

ax.set_xlabel("Rate estimate type")
ax.set_ylabel("Mutation rate (per bp, per gen.)")

handles, labels = ax.get_legend_handles_labels()

# When creating the legend, only use the first two elements
# to effectively remove the last two.
ax.legend(handles[3:], labels[3:], frameon=False, title="Haplotype at QTL")

sns.despine(ax=ax, top=True, right=True)

f.savefig(args.out, bbox_inches='tight')
