import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import pandas as pd
import argparse
import numpy as np

p = argparse.ArgumentParser()
p.add_argument("--mgp_msa")
args = p.parse_args()

msa = pd.read_csv(args.mgp_msa, sep='\t', header=None)

# remove strains that don't have associated singleton
# data from Dumont (2019)
not_in_dumont = ["CAST_EiJ", "LEWES_EiJ", "MOLF_EiJ", "PWK_PhJ", 
                 "SPRET_EiJ", "WSB_EiJ", "ZALENDE_EiJ"]

msa['reformatted_id'] = msa[0].apply(lambda s: s.lstrip('>').split('.')[0])

msa = msa[~msa['reformatted_id'].isin(not_in_dumont)]

f, ax = plt.subplots(figsize=(8, 2))

# strains that only have the "outside" three Mutyh mutations
withvar_some = ['BALB_cJ', 'BTBR_T+_Itpr3tf_J', 'BUB_BnJ', 'C3H_HeH', 'C3H_HeJ',
                'CBA_J', 'FVB_NJ', 'NOD_ShiLtJ', 'RF_J']

# strains that have all of the Mutyh mutations
withvar_all = ['A_J', 'DBA_1J', 'DBA_2J', 'ST_bJ']

# strains that have none of the Mutyh mutations
novar = ['C57BL_10J', 'C57BL_6NJ', 'C57BR_cdJ', 'C57L_J', 'C58_J', 'KK_HiJ',
         'NZB_B1NJ', 'NZO_HlLtJ', 'NZW_LacJ', 'SEA_GnJ']

# position of the Mutyh gene sequence start in mm10
mutyh_start = 116807723

# locations of five Mutyh missense mutations
mutyh_missense = [116814337, 116814393, 116815657, 116817415, 116817418]

idx2use = np.array([s - mutyh_start for s in mutyh_missense])

# get array of the nucleotides in the MSA
nucs = msa.values[:,1:-1].astype(np.int8)

samples = msa['reformatted_id']

# make a list of samples, ordered by their mutation status
ordered_smp = withvar_all
ordered_smp.extend(withvar_some)
ordered_smp.extend(novar)

smp2idx = dict(zip(samples, np.arange(len(samples))))
smp_idxs = np.array([smp2idx[s] for s in ordered_smp if s in smp2idx])

# define a custom colormap
myColors = ('gainsboro', 'cornflowerblue')
cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))

# make a heatmap object
ax = sns.heatmap(nucs[smp_idxs][:,idx2use].T, cmap=cmap, 
                    linewidths=.5, linecolor='k', ax=ax, cbar=True)

# Manually specify colorbar labelling after it's been generated
colorbar = ax.collections[0].colorbar
colorbar.set_ticks([0.5, 1.5])
colorbar.set_ticklabels(['REF', 'ALT'])

ax.set_xticks(np.arange(len(ordered_smp)) + 0.5)
ax.set_xticklabels([s.lstrip('>').split('.')[0].replace('_', '/').split('+')[0] \
                    for s in ordered_smp], fontsize=10, rotation=90)

def to_int(nuc):
    d = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
    return d[nuc] 
    
ax.set_yticklabels(['p.Gln5Arg', 'p.Arg24Cys', 'p.Ser69Arg', 
                    'p.Thr312Pro', 'p.Ser313Pro'], rotation=0, fontsize=10)
    
f.savefig('plots/mgp_wild_heatmap.eps', bbox_inches='tight')
