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
from collections import Counter

def clr(X):
    """
    perform a centered log-ratio transform
    """
    # the geometric mean acts as the center of the composition
    geom_mean = np.power(np.prod(X,axis=1),1/X.shape[1])
    return np.log(X / geom_mean[:,None])


p = argparse.ArgumentParser()
p.add_argument("--tidy_spectra", required=True, 
                    help="""tidy dataframe containing BXD mutation spectra""")
p.add_argument("--dumont_xls", required=True,
                    help="""Supplementary data from Dumont 2019""")
p.add_argument("--out", required=True,
                    help="""name of output plot""")
args = p.parse_args()
plt.rc('font', size=18)

dumont = pd.read_excel(args.dumont_xls, 
                       sheet_name="TableS3",
                       header=2)

dumont_filtered = dumont.query('Fraction_genome_in_IBD_Block >= 0.25')

## subset the dataframe to only include relevant info
## namely, the number of callable base pairs in each strain, dichotomized
## by the nucleotide at that base pair
dumont_filtered = dumont_filtered[['Strain', 
                          'C>A.2', 'C>G.2', 'C>T.2', 'T>A.2', 'T>C.2', 'T>G.2']]


## rename columns
colnames = ['strain', 'C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
dumont_filtered.columns = colnames

# change mutation types to batch the complement-collapsed
# types in the BXD dataframe
dumont_filtered.rename(columns={'T>A':'A>T', 'T>C':'A>G', 'T>G':'A>C'}, inplace=True)

# strains that only have the "outside" three Mutyh mutations
withvar_some = ['BALB_cJ', 'BTBR T+ Itpr3tf_J', 'BUB_BnJ', 'C3H_HeH', 'C3H_HeJ',
                'CBA_J', 'FVB_NJ', 'NOD_ShiLtJ', 'RF_J']

# strains that have all of the Mutyh mutations
withvar_all = ['A_J', 'DBA_1J', 'DBA_2J', 'ST_bJ']

# strains that have none of the Mutyh mutations
novar = ['C57BL_10J', 'C57BL_6NJ', 'C57BR_cdJ', 'C57L_J', 'C58_J', 'KK_HiJ',
         'NZB_B1NJ', 'NZO_HILtJ', 'NZW_LacJ', 'SEA_GnJ']

def dumont_strain_to_group(strain):
    if strain in withvar_all: return 'DBA_2J_like'
    elif strain in novar: return 'C57BL_6NJ_like'
    else: return 'intermediate'

## format strain names to be more readable
dumont_filtered['strain'] = dumont_filtered['strain'].apply(lambda x: x.replace('/', '_'))
dumont_filtered['sample_name'] = dumont_filtered['strain']

dumont_filtered['strain'] = dumont_filtered['strain'].apply(lambda x: dumont_strain_to_group(x))

# read in BXD mutation spectra
tidy_spectra = pd.read_csv(args.tidy_spectra)
tidy_fracs = tidy_spectra.query('estimate_type == "fraction"')
tidy_fracs = tidy_fracs[['bxd_strain_conv', 'haplotype_at_qtl', 'base_mut', 'estimate']]

# convert tidy to wide dataframe
wide_fracs = tidy_fracs.pivot(index=['bxd_strain_conv', 'haplotype_at_qtl'], 
                                columns='base_mut', values='estimate').reset_index()

# make a new column in the dataframe containing the sum of
# both types of C>T mutation fractions
wide_fracs['C>T.1'] = wide_fracs['C>T'] + wide_fracs['CpG>TpG']
wide_fracs.drop(columns=['C>T'], inplace=True)
wide_fracs.rename(columns={'C>T.1':'C>T', 'bxd_strain_conv': 'sample_name'}, inplace=True)

# change the strain column in the BXD data to just
# be the haplotype of the strain at the QTL
wide_fracs['strain'] = wide_fracs['haplotype_at_qtl']

muts = ['C>A', 'C>T', 'C>G', 'A>T', 'A>C', 'A>G']
new_colnames = muts
new_colnames.extend(['strain', 'sample_name'])
wide_fracs = wide_fracs[new_colnames]

# combine the dumont and BXD data into a single dataframe
combined = pd.concat([wide_fracs, dumont_filtered]).reset_index()
combined = combined.fillna(value=0)

combined = combined[combined['strain'].isin(['D', 'B', 'DBA_2J_like', 'C57BL_6NJ_like', 'intermediate'])]

muts = ['C>A', 'C>T', 'C>G', 'A>T', 'A>C', 'A>G']

X = combined[muts].values
y = combined['strain'].values

# get rid of any strains that have no mutation data
row_zero = np.unique(np.where(X == 0)[0])

X = np.delete(X, row_zero, axis=0)
y = np.delete(y, row_zero)

# perform centered log ratio transform on fraction data
#X = clr(X)

# perform PCA on data
#scaler = StandardScaler()
#X = scaler.fit_transform(X)
pca = PCA()
X_new = pca.fit_transform(X)

strains = set(y)
strain2color = {"D": "slateblue",
                "DBA_2J_like": "slateblue",
                "B": "lightgreen",
                "C57BL_6NJ_like": "lightgreen",
                "intermediate": "cornflowerblue"}

labels = [r'C$\to$A', r'C$\to$T', r'C$\to$G', r'A$\to$T', r'A$\to$C', r'A$\to$G']

sns.set_style('ticks')

f = plt.figure(figsize=(9, 7))
gs = gridspec.GridSpec(2, 1, height_ratios=[3, 2])
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])

x_vals = X_new[:,0]
y_vals = X_new[:,1]

colors = [strain2color[s] for s in y]

for idx in np.arange(x_vals.shape[0]):
    s = combined['sample_name'].values[idx]
    if combined['strain'].values[idx] in ("B", "D"): continue
    if x_vals[idx] > 0: print (s)

# plot each BXD and founder on PC plot

for strain_name in ["D", "DBA_2J_like", "B", "C57BL_6NJ_like", "intermediate"]:
    idxs = np.where(y == strain_name)[0]
    if strain_name in ["D", "B"]: marker = 'x'
    else: marker = 'o'
    ax1.scatter(x_vals[idxs], y_vals[idxs], color=strain2color[strain_name], 
        edgecolor='k', s=200, lw=2, marker=marker, label=strain_name)

exp_var = pca.explained_variance_ratio_

# plot PCA loadings
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

#legend_elements = [Patch(facecolor=c, edgecolor='k',
#                         label=m) for m,c in strain2color.items()]

#legend_elements = [Line2D([0], [0], color='slateblue', marker='x', label='D2 haplotype'),
#                   Line2D([0], [0], color='slateblue', marker='o', label='DBA/2J-like'),
#                   Line2D([0], [0], color='lightgreen', marker='x', label='B6 haplotype'),
#                   Line2D([0], [0], color='lightgreen', marker='o', label='C57BL/6NJ-like')]

ax2.set_ylim(ymin + (ymin / 2), ymax + (ymax / 2))
ax2.set_xlim(xmin + (xmin / 2), xmax + (xmax / 2))
ax2.set_xlabel("PC1 ({}%)".format(round(100 * exp_var[0]), 4), fontsize=18)
ax2.set_ylabel("PC2 ({}%)".format(round(100 * exp_var[1]), 3), fontsize=18)

#ax1.legend(handles=legend_elements, fontsize=12)
ax1.legend(fontsize=12)

f.tight_layout()

f.savefig(args.out, bbox_inches='tight')
