import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import seaborn as sns

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

# limit dataframe to only C>A mutations
df_ca = tidy_rates.query('base_mut == "C>A"')

# discretize the inbreeding time of each strain to be
# less than or greater than 20 generations
df_ca['Haplotype at QTL'] = df_ca['haplotype_at_qtl']
df_ca['Inbreeding time'] = df_ca['n_inbreeding_gens'].apply(lambda g: "< 20 generations" \
        if g < 20 else ">= 20 generations")

sns.set_style('ticks')

f, ax = plt.subplots()

sns.boxplot(x="Inbreeding time", y="estimate", hue="Haplotype at QTL", 
                        data=df_ca, ax=ax, color='w', fliersize=0)
sns.stripplot(x="Inbreeding time", y="estimate", hue="Haplotype at QTL", 
                        data=df_ca, ax=ax, palette='colorblind', dodge=True, 
                        edgecolor='k', lw=1.5)

ax.set_xlabel("Length of inbreeding in strain")
ax.set_ylabel(r'C$\to$A' + " mutation rate (per bp, per gen.)")

handles, labels = ax.get_legend_handles_labels()

# When creating the legend, only use the first two elements
# to effectively remove the last two.
ax.legend(handles[3:], labels[3:], frameon=False, title="Haplotype at QTL")

sns.despine(ax=ax, top=True, right=True)

f.savefig(args.out, bbox_inches='tight')
