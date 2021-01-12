import pandas as pd
import argparse
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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

tidy_rates.rename(columns={"epoch":"Epoch"}, inplace=True)

# ---
# generate figure 1a
# ---
sns.set_style('ticks')

f_1a, ax_1a = plt.subplots()

tidy_rates['muts_per_bp'] = tidy_rates['total_muts'] / tidy_rates['n_callable_bp']

sns.regplot(x="n_inbreeding_gens", y="muts_per_bp", data=tidy_rates, scatter=False,
        color='lightgrey', ax=ax_1a, line_kws={'zorder':0})
sns.scatterplot(x="n_inbreeding_gens", y="muts_per_bp", hue="Epoch", 
                palette="colorblind", edgecolor='w', data=tidy_rates, ax=ax_1a)


ax_1a.set_ylabel("Homozygous singleton mutations\nper callable base pair")
ax_1a.set_xlabel("Number of generations of inbreeding")

#ax_1a.legend(handles=legend_elements, frameon=False, title="Epoch of origin")

sns.despine(ax=ax_1a, right=True, top=True)
f_1a.savefig(args.outdir + 'main_paper/figure_1a.eps', bbox_inches='tight')

# ---
# generate figure 1b
# ---
f_1b, ax_1b = plt.subplots()


sns.boxplot(x="base_mut", y="estimate", linewidth=1, 
                color="w", fliersize=0, 
                data=tidy_spectra.query('estimate_type == "fraction" & n_inbreeding_gens >= 20'), ax=ax_1b)
sns.stripplot(x="base_mut", y="estimate", linewidth=1, 
                palette="colorblind", edgecolor='w', jitter=0.25,
                data=tidy_spectra.query('estimate_type == "fraction" & n_inbreeding_gens >= 20'), ax=ax_1b)

ax_1b.set_xticklabels(pd.unique(tidy_spectra['base_mut']))
ax_1b.set_ylabel("Fraction of homozygous singletons")
ax_1b.set_xlabel("Mutation type")

sns.despine(ax=ax_1b, right=True, top=True)
f_1b.savefig(args.outdir + 'main_paper/figure_1b.eps', bbox_inches='tight')
