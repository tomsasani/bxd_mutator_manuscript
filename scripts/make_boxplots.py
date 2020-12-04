import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns

plt.rc('font', size=12)

p = argparse.ArgumentParser()
p.add_argument("--tidy_rates")
p.add_argument("--tidy_fracs")
p.add_argument("--outdir")
args = p.parse_args()

mutation_rates_tidy = pd.read_csv(args.tidy_rates)
mutation_fracs_tidy = pd.read_csv(args.tidy_fracs)

sns.set_style('ticks')

f_1a, ax_1a = plt.subplots()

sns.boxplot(x="epoch", y="total_mutations", color='white', linewidth=1.5, 
                fliersize=0, data=mutation_rates_tidy, ax=ax_1a)
sns.swarmplot(x="epoch", y="total_mutations", palette='colorblind', linewidth=1, 
                edgecolor='k', data=mutation_rates_tidy, ax=ax_1a)

ax_1a.set_xticklabels(list(map(str, pd.unique(mutation_rates_tidy['epoch']))))
ax_1a.set_ylabel("Number of homozygous singletons")
ax_1a.set_xlabel("Epoch of origin")

sns.despine(ax=ax_1a, right=True, top=True)
f_1a.savefig(args.outdir + 'figure_1a.eps', bbox_inches='tight')

sns.set_style('ticks')

f_1b, ax_1b = plt.subplots()

sns.boxplot(x="epoch", y="rate", color='white', linewidth=1.5, 
                fliersize=0, data=mutation_rates_tidy, ax=ax_1b)
sns.swarmplot(x="epoch", y="rate", palette='colorblind', linewidth=1, 
                edgecolor='k', data=mutation_rates_tidy, ax=ax_1b)

ax_1a.set_xticklabels(list(map(str, pd.unique(mutation_rates_tidy['epoch']))))
ax_1b.set_ylabel("Mutation rate (per bp, per gen.)")
ax_1b.set_xlabel("Epoch of origin")

sns.despine(ax=ax_1b, right=True, top=True)
f_1b.savefig(args.outdir + '/figure_1b.eps', bbox_inches='tight')

sns.set_style('ticks')

f_1d, ax_1d = plt.subplots()

mutation_fracs_tidy = mutation_fracs_tidy.query('total_mutations >= 50')

flierprops = dict(marker='o', markersize=4, markeredgewidth=0, markerfacecolor='lightgrey')

sns.boxplot(x="base_mut", y="fraction", hue="epoch", linewidth=1.5, 
                palette="colorblind", flierprops=flierprops, data=mutation_fracs_tidy, ax=ax_1d)

ax_1d.set_xticklabels(pd.unique(mutation_fracs_tidy['base_mut']))
ax_1d.set_ylabel("Fraction of homozygous singletons")
ax_1d.set_xlabel("Mutation type")

sns.despine(ax=ax_1d, right=True, top=True)
f_1d.savefig(args.outdir + 'figure_1d.eps', bbox_inches='tight')
