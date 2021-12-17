from figure_gen_utils import *
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import argparse

p = argparse.ArgumentParser()
p.add_argument(
    "--tidy_spectra",
    required=True,
    help="""tidy dataframe containing BXD mutation spectra""",
)
p.add_argument(
    "--strain_metadata",
    required=True,
    help="""metadata about strains in Excel format.""",
)
args = p.parse_args()

tidy_spectra = pd.read_csv(args.tidy_spectra)
tidy_fracs = tidy_spectra.query(
    'estimate_type == "fraction" & base_mut == "C>A"')

# read in strain metadata and construct relevant dictionaries that
# map sample names to metadata
summary = pd.read_excel(args.strain_metadata)

# convert verbose filial generations to integer numbers of generations in the file
summary['gen_at_seq'] = summary['Generation at sequencing'].apply(
    get_generation)

summary['years_per_gen'] = summary.apply(calculate_years_per_gen, axis=1)

f, ax = plt.subplots()

sns.boxplot(
    x="Epoch",
    y="years_per_gen",
    data=summary.query("years_per_gen != -1"),
    ax=ax,
    color='white',
    fliersize=0,
)
sns.stripplot(
    x="Epoch",
    y="years_per_gen",
    data=summary.query("years_per_gen != -1"),
    ax=ax,
    palette="colorblind",
)

sns.despine(trim=True, ax=ax)
ax.set_ylabel("Generation time (years)")
f.savefig("plots/supp_figure_2a.eps")

import statsmodels.api as sm
from statsmodels.formula.api import ols

# Ordinary Least Squares (OLS) model
model = ols('years_per_gen ~ C(Epoch)',
            data=summary.query("years_per_gen != -1")).fit()
anova_table = sm.stats.anova_lm(model, typ=2)

print(np.mean(summary.query("years_per_gen != -1").years_per_gen))
print(anova_table)

strain2gentime = dict(zip(summary['bam_name'], summary['years_per_gen']))
tidy_fracs['years_per_gen'] = tidy_fracs['bxd_strain'].apply(
    lambda b: strain2gentime[b])
model = ols('estimate ~ years_per_gen + C(epoch)',
            data=tidy_fracs.query("years_per_gen != -1")).fit()
anova_table = sm.stats.anova_lm(model, typ=2)
print(anova_table)

f, ax = plt.subplots()

sns.regplot(
    x="years_per_gen",
    y="estimate",
    color="cornflowerblue",
    data=tidy_fracs,
    scatter_kws={"edgecolor": "k"},
)

ax.set_xlim(0.15, 0.75)

ax.set_xlabel("Generation time (years)")
ax.set_ylabel("C" + r"$\rightarrow$" + "A singleton fraction")

f.savefig("plots/supp_figure_2b.eps")
