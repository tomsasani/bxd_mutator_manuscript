import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import seaborn as sns

plt.rc('font', size=12)

p = argparse.ArgumentParser()
p.add_argument("--tidy_rates", required=True, 
                    help="""file containing adjusted mutation rates for \
                            each sample in tidy, long-form format""")
p.add_argument("--tidy_spectra", required=True,
                    help="""file containing mutation spectra for \
                            each sample in tidy, long-form format""")
p.add_argument("--outdir", required=True,
                    help="""name of directory to store output figures""")
args = p.parse_args()

# read in tidy data
tidy_rates = pd.read_csv(args.tidy_rates)
tidy_spectra = pd.read_csv(args.tidy_spectra)

tidy_rates = tidy_rates.astype({'epoch': 'int'})
tidy_spectra = tidy_spectra.astype({'epoch': 'int'})

print ('total of {} mutations in {} strains'.format(np.sum(tidy_rates['total_muts']), len(pd.unique(tidy_rates['bxd_strain_conv']))))

# ---
# generate figure 1a
# ---
sns.set_style('ticks')

f_1a, ax_1a = plt.subplots()

sns.boxplot(x="epoch", y="total_muts", color='white', linewidth=1.5, 
                fliersize=0, data=tidy_rates, ax=ax_1a)
sns.stripplot(x="epoch", y="total_muts", palette='colorblind', linewidth=1, 
                edgecolor='k', data=tidy_rates, ax=ax_1a)

ax_1a.set_xticklabels(list(map(str, pd.unique(tidy_rates['epoch']))))
ax_1a.set_ylabel("Number of homozygous singletons")
ax_1a.set_xlabel("Epoch of origin")

sns.despine(ax=ax_1a, right=True, top=True)
f_1a.savefig(args.outdir + 'figure_1a.eps', bbox_inches='tight')

# ---
# generate figure 1b
# ---
sns.set_style('ticks')

f_1b, ax_1b = plt.subplots()

sns.boxplot(x="epoch", y="rate", color='white', linewidth=1.5, 
                fliersize=0, data=tidy_rates, ax=ax_1b)
sns.stripplot(x="epoch", y="rate", palette='colorblind', linewidth=1, 
                edgecolor='k', data=tidy_rates, ax=ax_1b)

ax_1a.set_xticklabels(list(map(str, pd.unique(tidy_rates['epoch']))))
ax_1b.set_ylabel("Mutation rate (per bp, per gen.)")
ax_1b.set_xlabel("Epoch of origin")

sns.despine(ax=ax_1b, right=True, top=True)
f_1b.savefig(args.outdir + '/figure_1b.eps', bbox_inches='tight')

sns.set_style('ticks')

# ---
# generate figure 1d
# ---
f_1d, ax_1d = plt.subplots()

#tidy_spectra['base_mut_new'] = tidy_spectra['base_mut'].apply(lambda m: m.replace(">", r'$\to$'))

flierprops = dict(marker='o', markersize=0, markeredgewidth=0, markerfacecolor='lightgrey')

print (tidy_spectra.query('estimate_type == "fraction" & base_mut == "CpG>TpG" & estimate > 0.3'))

sns.boxplot(x="base_mut", y="estimate", linewidth=1, 
                color='w', fliersize=0, #=flierprops, 
                data=tidy_spectra.query('estimate_type == "fraction"'), ax=ax_1d)
sns.stripplot(x="base_mut", y="estimate", linewidth=1, 
                palette='colorblind', #flierprops=flierprops, 
                data=tidy_spectra.query('estimate_type == "fraction"'), ax=ax_1d)
# perform ANOVA to test for differences in spectra across epochs
import statsmodels.api as sm
from statsmodels.formula.api import ols

tidy_fractions = tidy_spectra.query('estimate_type == "fraction"')
model = ols('estimate ~ C(base_mut) + C(epoch)', data=tidy_fractions).fit()
print (sm.stats.anova_lm(model, typ=2))

ax_1d.set_xticklabels(pd.unique(tidy_spectra['base_mut']))
ax_1d.set_ylabel("Fraction of homozygous singletons")
ax_1d.set_xlabel("Mutation type")

sns.despine(ax=ax_1d, right=True, top=True)
f_1d.savefig(args.outdir + 'figure_1d.eps', bbox_inches='tight')
