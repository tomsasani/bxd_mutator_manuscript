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
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib import gridspec
from collections import Counter

plt.rc('font', size=18)

p = argparse.ArgumentParser()
p.add_argument("--wild_singleton_vars", required=True, nargs="*",
                    help="""paths to per-chromosome CSVs of wild mouse singletons""")
p.add_argument("--out", required=True,
                    help="""name of output plot""")
args = p.parse_args()

# combine the per-chromosome csvs of mutation counts
df = None
for fh in args.wild_singleton_vars:
    if df is None:
        df = pd.read_csv(fh).set_index("Unnamed: 0")
    else:
        df = df.add(pd.read_csv(fh).set_index("Unnamed: 0"))
df = df.reset_index()

df['species'] = df['Unnamed: 0'].apply(lambda s: s.split('_')[0])

# map mutations to indices
muts = list(df)[1:-1]
mut2idx = dict(zip(muts, range(len(muts))))

colnames = muts.copy()
colnames.append('species')

df = df[colnames]

X = df[muts].values
y = df['species'].values

# convert mutation counts to fractions
X_sums = np.sum(X, axis=1)
X = X / X_sums[:,None]

# perform PCA on data
scaler = StandardScaler()
#X = scaler.fit_transform(X)
pca = PCA()
X_new = pca.fit_transform(X)

sns.set_style('ticks')

# create plot object
f, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

# get PC1 and PC2
x_vals = X_new[:,0]
y_vals = X_new[:,1]

colors = sns.color_palette('colorblind', len(set(y)))

labels = [r'A$\to$T', r'A$\to$C', r'A$\to$G', r'C$\to$A', r'C$\to$T', r'C$\to$G']

for cat_i,cat in enumerate(["Mmd", "Mmc", "Mmm", "Ms"]):
    idxs = np.where(y == cat)[0]

    # plot samples
    ax1.scatter(x_vals[idxs], y_vals[idxs], color=colors[cat_i], 
            edgecolor='k', s=200, lw=2, label=cat)

# plot PCA loadings
coef = np.transpose(pca.components_[0:2, :])
ymin, ymax = 1, -1
xmin, xmax = 1, -1
for i in range(coef.shape[0]):
    ax2.arrow(0, 0, coef[i,0], coef[i,1], color = 'k', lw=2)
    t_x, t_y = coef[i,0] * 1.15, coef[i,1] * 1.15
    ax2.text(t_x, t_y, labels[i], color = 'k', ha = 'center', va = 'center')


ax1.legend()

exp_var = pca.explained_variance_ratio_

ax1.set_xlabel("PC1 ({}%)".format(round(100 * exp_var[0]), 4), fontsize=18)
ax1.set_ylabel("PC2 ({}%)".format(round(100 * exp_var[1]), 3), fontsize=18)

f.tight_layout()

f.savefig(args.out, bbox_inches='tight')