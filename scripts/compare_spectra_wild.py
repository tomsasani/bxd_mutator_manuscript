import argparse
import scipy.stats as ss
import itertools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from utils import revcomp, combine_chr_df
from mutation_comparison import mutation_comparison
import glob

p = argparse.ArgumentParser()
p.add_argument("--wild_singleton_vars", required=True, nargs="*",
                    help="""paths to per-chromosome CSVs of wild mouse singletons""")
p.add_argument("--out", required=True,
                    help="""name of output plot""")
args = p.parse_args()

# combine singleton counts from all chromosomes
df = None
for fh in glob.glob(args.wild_singleton_vars):
    if df is None:
        df = pd.read_csv(fh).set_index("Unnamed: 0")
    else:
        df = df.add(pd.read_csv(fh).set_index("Unnamed: 0"))

# add columns indicating species or subspecies of each strain
df = df.reset_index()
df['species'] = df['Unnamed: 0'].apply(lambda s: s.split('_')[0])
df['subspecies'] = df['Unnamed: 0'].apply(lambda s: s.split('_')[1])

#df = df.query('subspecies != "IRA"')

muts = list(df)[1:-1]

mut2idx = dict(zip(muts, range(len(muts))))

# for every comparison of domesticus to the other wild strains,
# generate heatmaps comparing mutation fractions
for combo in (["Mmd", "Mmm"], ["Mmd", "Mmc"], ["Mmd", "Ms"]):

    a_lab, b_lab = combo

    if a_lab == b_lab: continue

    a = df[df['species'] == a_lab]
    b = df[df['species'] == b_lab]

    subset_0 = np.sum(a.values[:,1:-2], axis=0)
    subset_1 = np.sum(b.values[:,1:-2], axis=0)

    mutation_comparison(subset_0, subset_1, mut2idx=mut2idx, 
                    outname="plots/mutation_enrichment_wild.{}_{}.eps".format(a_lab, b_lab), 
                    nmer4norm=None, plot_type="heatmap",
                    title=None)