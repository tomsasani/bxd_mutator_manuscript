import pandas as pd
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

# ---
# generate figure 1a
# ---
sns.set_style('ticks')

f_1a, ax_1a = plt.subplots()

sns.boxplot(x="epoch", y="total_muts", color='white', linewidth=1.5, 
                fliersize=0, data=tidy_rates, ax=ax_1a)
sns.swarmplot(x="epoch", y="total_muts", palette='colorblind', linewidth=1, 
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
sns.swarmplot(x="epoch", y="rate", palette='colorblind', linewidth=1, 
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

flierprops = dict(marker='o', markersize=4, markeredgewidth=0, markerfacecolor='lightgrey')

sns.boxplot(x="base_mut", y="estimate", hue="epoch", linewidth=1, 
                palette="colorblind", flierprops=flierprops, 
                data=tidy_spectra.query('estimate_type == "fraction"'), ax=ax_1d)

ax_1d.set_xticklabels(pd.unique(tidy_spectra['base_mut']))
ax_1d.set_ylabel("Fraction of homozygous singletons")
ax_1d.set_xlabel("Mutation type")

ax_1d.legend(frameon=False, title="Epoch of origin", loc="upper left", fontsize=10)

sns.despine(ax=ax_1d, right=True, top=True)
f_1d.savefig(args.outdir + 'figure_1d.eps', bbox_inches='tight')
