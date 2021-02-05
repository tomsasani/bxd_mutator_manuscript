import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns 
import argparse

p = argparse.ArgumentParser()
p.add_argument("--strain_vars")
p.add_argument("--out")
args = p.parse_args()

df = pd.read_csv(args.strain_vars, 
                    names=["interval", "strain", "gt", "ref", "alt"])

# remove strains that don't have associated singleton
# data from Dumont (2019)
not_in_dumont = ["CAST_EiJ", "LEWES_EiJ", "MOLF_EiJ", "PWK_PhJ", 
                 "SPRET_EiJ", "WSB_EiJ", "ZALENDE_EiJ"]

df = df[~df['strain'].isin(not_in_dumont)]

samps = pd.unique(df['strain'])
intervals = pd.unique(df['interval'])

out_a = np.zeros((len(intervals), len(samps)))

smp2idx = dict(zip(samps, range(len(samps))))
int2idx = dict(zip(intervals, range(len(intervals))))

for i,row in df.iterrows():
    i_i = int2idx[row['interval']]
    strain = row['strain']
    
    s_i = smp2idx[row['strain']]
    out_a[i_i, s_i] = row['gt']

sns.set_style('ticks')
f, ax = plt.subplots(figsize=(16,4))

cmap = ['gainsboro', 'royalblue']

interval2mut = dict(zip(intervals, ["p.Gln5Arg",
                                   "p.Arg24Cys",
                                   "p.Ser69Arg",
                                   "p.Thr312Pro",
                                   "p.Thr313Pro"]))

# get the number of Mutyh mutations in each strain, and sort
# the output heatmap in that order
n_mutations = np.sum(out_a, axis=0)
sorted_idxs = np.argsort(n_mutations)[::-1]

out_a = out_a[:,sorted_idxs]

sns.heatmap(out_a, ax=ax, cmap=cmap, linecolor='k', linewidth=1)

colorbar = ax.collections[0].colorbar
colorbar.set_ticks([0.5, 1.5,])
colorbar.set_ticklabels(['C57BL/6J allele', 'DBA/2J allele'])

samps_reform = np.array([s.replace('_', '/') for s in samps])
ax.set_xticks(np.arange(len(samps_reform)) + 0.5)
ax.set_xticklabels(samps_reform[sorted_idxs], rotation=90, fontsize=16)

ax.set_yticks(np.arange(len(intervals)) + 0.5)
ax.set_yticklabels([interval2mut[it] for it in intervals], fontsize=16, rotation=0)

f.tight_layout()
f.savefig(args.out, bbox_inches='tight')

