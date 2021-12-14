from typing import ValuesView
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import argparse
import numpy as np

p = argparse.ArgumentParser()
p.add_argument(
    "--wild_sfs",
    nargs="*",
    required=True,
    help="""paths to wild SFS files""",
)
p.add_argument(
    "--out",
    required=True,
    help="""name of output file""",
)
args = p.parse_args()

species = [a.split('.')[-2] for a in args.wild_sfs]

f, ax = plt.subplots()

colors = sns.color_palette('colorblind', len(args.wild_sfs))

for si, species_name in enumerate(species):
    out_a = None
    mut2idx = None
    af = None
    ac = None
    for sfs_i, sfs in enumerate(args.wild_sfs):
        if sfs.split('.')[-2] != species_name: continue
        df = pd.read_csv(sfs, sep='\t')
        if mut2idx is None:
            mut2idx = dict(zip(list(df)[1:], range(len(list(df)[1:]))))
            ac = df.values[:, 0]
            af = ac / np.max(ac)
        mut_counts = df.values[:, 1:]
        if out_a is None:
            out_a = mut_counts
        else:
            out_a += mut_counts
    mut_sums = np.sum(out_a, axis=1)
    mut_fracs = out_a / mut_sums[:, None]
    ax.plot(
        af,
        mut_fracs[:, mut2idx["C>A"]],
        label=species_name,
        color=colors[si],
        lw=2,
    )

ax.set_xlabel("Allele frequency")
ax.set_ylabel("C>A mutation fraction")
ax.legend()

f.tight_layout()
f.savefig(args.out)