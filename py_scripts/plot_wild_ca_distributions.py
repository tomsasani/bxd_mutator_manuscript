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

plt.rcParams.update({'font.size': 16})

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
df['cast_dom'] = df['species'].apply(lambda s: "Castaneus/Domesticus" if s in ["Mmc", "Mmd"] else "Musculus/Spretus")

header = list(df)

# make figure object and define colors
f, axarr = plt.subplots(len(df.kmer.unique()), figsize=(14, 20), sharey=True)
colors = sns.color_palette('colorblind', 4)

windows = df.region.unique()
# loop over all of the mutations 
for i, mut in enumerate(df.kmer.unique()):
    xs, ys = defaultdict(list), defaultdict(list)
    gts = []
    # loop over all of the windows
    for window, window_df in df.groupby("region"):
        # figure out the bp position
        pos = int(float(window.split(':')[1].split('-')[0]))

        conting = []
        # get the total counts of the focal mutation and the other
        # mutation types in each "genotype class"
        for genotype, genotype_df in window_df.groupby('cast_dom'):
            fore = np.sum(genotype_df.query("kmer == @mut")["count"])
            back = np.sum(genotype_df.query("kmer != @mut")["count"])
            # add each genotype's "foreground" and "background" counts
            # to a contingency table
            conting.append([fore, back])
            if genotype not in gts:
                gts.append(genotype)
        # ignore windows where one of the 
        if any([0 in el for el in conting]): continue
        # do the chi2 contingency test to figure out if the focal mutation
        # type is enriched in one genotype class vs the other
        _, p, _, _ = ss.chi2_contingency(conting)

        if p < 0.05: print (mut, window, p)

        # loop over the genotypes and add to the x/y value 
        # dictionaries for each
        for gi, g in enumerate(gts):
            xs[g].append(pos)
            ys[g].append(conting[gi][0] / np.sum(conting[gi]))

    
    for gi, g in enumerate(xs):
        #axarr[i].scatter(xs[g], ys[g], color=colors[gi], ec='k', label=g)
        axarr[i].plot(xs[g], ys[g], color=colors[gi], label=g, lw=2)

    axarr[i].set_title(mut)
    axarr[i].legend()
    axarr[i].set_xlabel("Position on chromosome 4")
    axarr[i].set_ylabel("Mutation fraction")

f.tight_layout()
f.savefig(args.out)
