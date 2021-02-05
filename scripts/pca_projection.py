import argparse
from matplotlib.lines import Line2D
import scipy.stats as ss
import itertools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import scale, normalize
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib import gridspec

p = argparse.ArgumentParser()
p.add_argument("--tidy_spectra")
p.add_argument("--dumont_xls")
p.add_argument("--out")
args = p.parse_args()
plt.rc('font', size=18)

dumont = pd.read_excel(args.dumont_xls, 
                       sheet_name="TableS3",
                       header=2)

## subset the dataframe to only include relevant info
## namely, the number of callable base pairs in each strain, dichotomized
## by the nucleotide at that base pair
dumont_filtered = dumont[['Strain', 
                          'C>A.2', 'C>G.2', 'C>T.2', 'T>A.2', 'T>C.2', 'T>G.2']]

## rename columns
colnames = ['strain', 'C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
dumont_filtered.columns = colnames

## format strain names to be more readable
dumont_filtered['strain'] = dumont_filtered['strain'].apply(lambda x: x.replace('/', '_'))

tidy_spectra = pd.read_csv(args.tidy_spectra)

tidy_fracs = tidy_spectra.query('estimate_type == "fraction"')

tidy_fracs = tidy_fracs.query('n_inbreeding_gens >= 20')
tidy_fracs = tidy_fracs[['bxd_strain_conv', 'haplotype_at_qtl', 'base_mut', 'estimate']]

wide_fracs = tidy_fracs.pivot(index=['bxd_strain_conv', 'haplotype_at_qtl'], columns='base_mut', values='estimate').reset_index()

wide_fracs['C>T.1'] = wide_fracs['C>T'] + wide_fracs['CpG>TpG']
wide_fracs.drop(columns=['C>T'], inplace=True)
wide_fracs.rename(columns={'C>T.1':'C>T'}, inplace=True)

wide_fracs['strain'] = wide_fracs['haplotype_at_qtl'].apply(lambda h: 'D2' if h == 1 else 'B6')

wide_fracs.rename(columns={'A>T':'T>A', 'A>G':'T>C', 'A>C':'T>G'}, inplace=True)

wide_fracs = wide_fracs[colnames]

combined = pd.concat([wide_fracs, dumont_filtered]).reset_index()

combined = combined.fillna(value=0)

combined = combined[combined['strain'].isin(['D2', 'B6', 'DBA_2J', 'C57BL_6NJ'])]

combined.to_csv("csv/pca.csv", index=False)

muts = ['C>A', 'C>T', 'C>G', 'T>A', 'T>G', 'T>C']
X = combined[muts].values
y = combined['strain'].values

def clr(X):
    # the geometric mean acts as the center of the composition
    geom_mean = np.power(np.prod(X,axis=1),1/X.shape[1])
    return np.log(X / geom_mean[:,None])

row_zero = np.unique(np.where(X == 0)[0])

print (row_zero)
print (X[row_zero])

X = np.delete(X, row_zero, axis=0)
y = np.delete(y, row_zero)

X = clr(X)

strains = set(y)
strain2color = {"D2": "slateblue",
                "DBA_2J": "slateblue",
                "B6": "lightgreen",
                "C57BL_6NJ": "lightgreen"}


pca = PCA()
X_new = pca.fit_transform(X)

labels = [r'C$\to$A', r'C$\to$T', r'C$\to$G', r'T$\to$A', r'T$\to$G', r'T$\to$C']

sns.set_style('ticks')

f = plt.figure(figsize=(10, 8))
gs = gridspec.GridSpec(2, 1, height_ratios=[3, 2])
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])#, sharex=ax1)

x_vals = X_new[:,0]
y_vals = X_new[:,1]

colors = [strain2color[s] for s in y]

ax1.scatter(x_vals[:-2], y_vals[:-2], c=colors[:-2], edgecolor='k', s=200, lw=2, marker='x')
ax1.scatter(x_vals[-2:], y_vals[-2:], c=colors[-2:], edgecolor='k', s=200, lw=2, marker='o')

legend_elements = [Patch(facecolor=c, edgecolor='k',
                         label=m) for m,c in strain2color.items()]

legend_elements = [Line2D([0], [0], color='slateblue', marker='x', label='D2 haplotype'),
                   Line2D([0], [0], color='slateblue', marker='o', label='DBA/2J'),
                   Line2D([0], [0], color='lightgreen', marker='x', label='B6 haplotype'),
                   Line2D([0], [0], color='lightgreen', marker='o', label='C57BL/6NJ')]

exp_var = pca.explained_variance_ratio_

coef = np.transpose(pca.components_[0:2, :])
ymin, ymax = 1, -1
xmin, xmax = 1, -1
for i in range(coef.shape[0]):
    ax2.arrow(0, 0, coef[i,0], coef[i,1], color = 'k', lw=2)
    t_x, t_y = coef[i,0] * 1.15, coef[i,1] * 1.15
    ax2.text(t_x, t_y, labels[i], color = 'k', ha = 'center', va = 'center')
    if t_y > ymax: ymax = t_y 
    if t_y < ymin: ymin = t_y 
    if t_x > xmax: xmax = t_x 
    if t_x < xmin: xmin = t_x 

ax2.set_ylim(ymin + (ymin / 2), ymax + (ymax / 2))
ax2.set_xlim(xmin + (xmin / 2), xmax + (xmax / 2))
ax2.set_xlabel("PC1 ({}%)".format(round(100 * exp_var[0]), 4), fontsize=18)
ax2.set_ylabel("PC2 ({}%)".format(round(100 * exp_var[1]), 3), fontsize=18)

ax1.legend(handles=legend_elements)

f.tight_layout()

f.savefig(args.out, bbox_inches='tight')
