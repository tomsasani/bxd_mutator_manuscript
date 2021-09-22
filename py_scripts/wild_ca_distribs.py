import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats as ss
import seaborn as sns
from collections import defaultdict
import argparse

p = argparse.ArgumentParser()
p.add_argument("--windowed_vars")
args = p.parse_args()

plt.rcParams.update({'font.size': 16})

#df = pd.read_csv(args.windowed_vars)
df = pd.read_csv("wild.singletons.csv")
groupby_cols = [
    'kmer',
    'region',
    'count',
    'frac',
    'sample',
]

df = df.drop_duplicates(groupby_cols)

print (df.head())

df['species'] = df['sample'].apply(lambda s: s.split('_')[0])

df['cast_dom'] = df['species'].apply(lambda s: "Castaneus/Domesticus" if s in ["Mmc", "Mmd"] else "Musculus/Spretus")

header = list(df)

df = df.query('genotype != -1')

f, axarr = plt.subplots(6, figsize=(14, 20), sharey=True)

windows = df.region.unique()

for i, mut in enumerate(df.kmer.unique()):

    gts = []

    idx = 0

    xs, ys = defaultdict(list), defaultdict(list)
    ps = []
    for window, window_df in df.groupby("region"):

        pos = int(float(window.split(':')[1].split('-')[0]))
        conting = []
        ca_fracs = []

        for genotype, genotype_df in window_df.groupby('cast_dom'):
            fore = np.sum(genotype_df.query("kmer == @mut")["count"])
            back = np.sum(genotype_df.query("kmer != @mut")["count"])

            conting.append([fore, back])
            if genotype not in gts:
                gts.append(genotype)

        if any([0 in el for el in conting]): continue
        #if any([np.sum(c) < 100 for c in conting]): continue

        _, p, _, _ = ss.chi2_contingency(conting)

        trans_p = np.log10(p) * -1
        print(conting, p)
        if trans_p > 3:
            ps.append("firebrick")
        else:
            ps.append("grey")

        colors = sns.color_palette('colorblind', len(gts))
        #axarr[i].scatter(idx, np.log10(p) * -1, color='grey', ec='k')
        for gi, g in enumerate(gts):
            xs[g].append(pos)
            ys[g].append(conting[gi][0] / np.sum(conting[gi]))
            #axarr[i].scatter(pos, conting[gi][0] / np.sum(conting[gi]), color=colors[gi], ec='k')
            #alpha=(np.log10(p) * -1) / 100. )#, label=g)
    for gi, g in enumerate(xs):
        axarr[i].scatter(xs[g], ys[g], color=colors[gi], ec='k', label=g)
    axarr[i].set_title(mut)
    axarr[i].legend()
    axarr[i].set_xlabel("Position on chromosome 4")
    axarr[i].set_ylabel("Mutation fraction")
#ax.set_xticks(np.arange(len(windows))[::5])
#ax.set_xticklabels([w.split(':')[-1].split('-')[0] for w in windows][::5], rotation=90)
f.tight_layout()
f.savefig('o.png')
