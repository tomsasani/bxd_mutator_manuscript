import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats as ss
import seaborn as sns
from collections import defaultdict
import argparse

p = argparse.ArgumentParser()
p.add_argument("--windowed_vars")
p.add_argument("--out")
args = p.parse_args()

plt.rcParams.update({'font.size': 20})

df = pd.read_csv(args.windowed_vars)

groupby_cols = [
    'kmer',
    'region',
    'count',
    'frac',
    'sample',
]

df = df.drop_duplicates(groupby_cols)

# add species and "genotype" columns
# the "genotype" column is essentially a binary, saying whether
# the strain is Cast/Dom or Mus/Spret
df['species'] = df['sample'].apply(lambda s: s.split('_')[0])

header = list(df)

# make figure object and define colors
#f, axarr = plt.subplots(len(df.kmer.unique()), figsize=(14, 20), sharey=True)
f, ax = plt.subplots(figsize=(12, 8))

colors = sns.color_palette('colorblind', 4)

out_df = []

windows = df.region.unique()
# loop over all of the mutations
for i, mut in enumerate(df.kmer.unique()):
    xs, ys = defaultdict(list), defaultdict(list)
    gts = []
    # loop over all of the windows
    for window, window_df in df.groupby("region"):
        # figure out the bp position
        pos = int(float(window.split(':')[1].split('-')[0]))

        if pos < 114800000 or pos > 118300000: continue

        conting = []
        # get the total counts of the focal mutation and the other
        # mutation types in each "genotype class"
        for genotype, genotype_df in window_df.groupby('species'):
            fore = np.sum(genotype_df.query("kmer == @mut")["count"])
            back = np.sum(genotype_df.query("kmer != @mut")["count"])
            # add each genotype's "foreground" and "background" counts
            # to a contingency table
            conting.append([fore, back])
            if genotype not in gts:
                gts.append(genotype)

        # ignore windows where one of the
        if any([0 in el for el in conting]): continue

        gt2idx = dict(zip(gts, range(len(gts))))

        cast_counts, dom_counts = conting[gt2idx['Mmc']], conting[gt2idx['Mmd']]
        spret_counts, musc_counts = conting[gt2idx['Ms']], conting[gt2idx['Mmm']]

        for mouse in gts:
            counts = conting[gt2idx[mouse]]
            frac = counts[0] / sum(counts)
            out_dict = {'window': window, 'mutation_type': mut, 'Species': mouse, 'fraction': frac}
            out_df.append(out_dict)

out_df = pd.DataFrame(out_df)
sns.stripplot(x="mutation_type", y="fraction", hue="Species", dodge=True, data=out_df, ax=ax, palette='colorblind')

import statsmodels.api as sm
from statsmodels.formula.api import ols
model = ols('fraction ~ C(mutation_type) + C(Species)', data=out_df).fit()
anova_table = sm.stats.anova_lm(model, typ=2)
print (anova_table)

for species, species_df in out_df.query('mutation_type == "C>A"').groupby('Species'):
    print (species, ss.median_abs_deviation(species_df.fraction.values))


ax.set_xlabel("Mutation type")
ax.set_ylabel("Singleton fractions in 50-kbp windows\nacross the chr4 QTL interval")

f.tight_layout()
f.savefig(args.out)