import argparse
import scipy.stats as ss
import itertools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

p = argparse.ArgumentParser()
p.add_argument("--dumont_xls", required=True,
                    help="""strain-private substitution data from Dumont et al. (2019) \
                        as it appears in Table 3 of the Supplementary Data""")
p.add_argument("--out", required=True,
                    help="""name of output plot""")
p.add_argument("-mutation_type", default="C>A")
args = p.parse_args()

dumont = pd.read_excel(args.dumont_xls, 
                       sheet_name="TableS3",
                       header=2)

# subset the dataframe to only include relevant info
# namely, the number of callable base pairs in each strain, dichotomized
# by the nucleotide at that base pair
dumont_filtered = dumont[['Strain', 'nA', 'nC', 'nG', 'nT',
                          'C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']]

# rename columns
dumont_filtered.columns = ['strain', 'nA', 'nC', 'nG', 'nT',
                           'C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']

# format strain names to be more readable
dumont_filtered['strain'] = dumont_filtered['strain'].apply(lambda x: x.replace('/', '_'))

# map mutations to corresponding indices
muts = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
mut2idx = dict(zip(muts, range(len(muts))))

# map samples to corresponding indices
smp2idx = dict(zip(dumont_filtered['strain'], 
                   range(len(dumont_filtered['strain']))))

# get the counts of singletons of each mutation type
# in each strain
counts = dumont_filtered.values[:,5:11]

# get the counts of callable base pairs of each nucleotide
# type in each strain
denoms = dumont_filtered.values[:,1:5]

# strains that only have the "outside" three Mutyh mutations
withvar_some = ['BALB_cJ', 'BTBR T+ Itpr3tf_J', 'BUB_BnJ', 'C3H_HeH', 'C3H_HeJ',
                'CBA_J', 'FVB_NJ', 'NOD_ShiLtJ', 'RF_J']

loner = ["I_LnJ"]

# strains that have all of the Mutyh mutations
withvar_all = ['A_J', 'DBA_1J', 'DBA_2J', 'ST_bJ']

# strains that have none of the Mutyh mutations
novar = ['C57BL_10J', 'C57BL_6NJ', 'C57BR_cdJ', 'C57L_J', 'C58_J', 'KK_HiJ',
         'NZB_B1NJ', 'NZO_HILtJ', 'NZW_LacJ', 'SEA_GnJ']

# map strains to the "category" of Mutyh mutations they belong to
cat2strain = {'DBA-like': withvar_all, 'intermediate': withvar_some, 
                'I_LnJ': loner, 'C57-like': novar}

f, axarr = plt.subplots(1,2, sharey=True, figsize=(14, 4))

sns.set_style('ticks')

# do every pairwise comparison of spectra in strains with each
# configuration of Mutyh mutations
for cat_i,cat in enumerate([("DBA-like", "intermediate"), ("intermediate", "C57-like")]):

    a, b = cat

    a_idx = np.array([smp2idx[s] for s in cat2strain[a] if s in smp2idx])
    b_idx = np.array([smp2idx[s] for s in cat2strain[b] if s in smp2idx])
    
    # subsets of mutation counts for strains in each
    # mutation "category"
    subset_0 = counts[a_idx, :]
    subset_1 = counts[b_idx, :]

    # get the callable number of basepairs for each 
    # nucleotide for the samples in either category
    subset_0_denom = denoms[a_idx, :]
    subset_1_denom = denoms[b_idx, :]

    # sum the mutation counts of all types for strains
    # in each category
    subset_0_totals = np.sum(subset_0, axis=0)
    subset_1_totals = np.sum(subset_1, axis=0)
  
    # get the total number of callable base pairs (across all
    # nucleotides) in either subset
    subset_0_denom_totals = np.sum(subset_0_denom, axis=0)
    subset_1_denom_totals = np.sum(subset_1_denom, axis=0)

    nuc2denom = {'A':0, 'C':1, 'G':2, 'T':3}

    nuc2comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

    ind = np.arange(len(muts))

    a_adj, b_adj = [], []

    # we want to perform a Chi-square test to compare mutation
    # spectra in subset 0 vs. subset 1, but we need to first
    # adjust the counts of each mutation type in each subset by
    # the number of A, C, T, or G base pairs that met filtering
    # criteria across all strains in either subset.
    for i,mut in enumerate(muts):
        base_nuc = mut.split('>')[0]

        # first, sum the total number of mutations of type `mut`
        # observed in each subset
        a_fore = subset_0_totals[i]
        b_fore = subset_1_totals[i]

        # then get thethe sum of callable base pairs that are
        # X nucleotides across strains in subset 0, and same for subset 1.
        # we add the sum of callable base pairs that
        # are the complement of X, as well
        a_denom = subset_0_denom_totals[nuc2denom[base_nuc]]
        a_denom += subset_0_denom_totals[nuc2denom[nuc2comp[base_nuc]]]
        b_denom = subset_1_denom_totals[nuc2denom[base_nuc]]
        b_denom += subset_1_denom_totals[nuc2denom[nuc2comp[base_nuc]]]

        # we now adjust the counts of each mutation type in each strain
        if a_denom > b_denom: 
            adj = b_denom / a_denom
            a_fore = a_fore * adj
        elif b_denom > a_denom:
            adj = a_denom / b_denom
            b_fore = b_fore * adj

        a_adj.append(a_fore)
        b_adj.append(b_fore)

    # for every mutation type, perform a chi-square test of
    # enrichment of that mutation type's counts in subset 0 vs. subset 1
    mut_i = mut2idx[args.mutation_type]
    a_fore = a_adj[mut_i]
    a_back = sum(a_adj) - a_fore
    b_fore = b_adj[mut_i]
    b_back = sum(b_adj) - b_fore
    
    a_adj_fracs = np.array(a_adj) / sum(a_adj)
    b_adj_fracs = np.array(b_adj) / sum(b_adj)
    
    ind = np.arange(len(mut2idx))

    axarr[cat_i].bar(ind - 0.2, a_adj_fracs, 0.4, 
                     label=cat[0], color='cornflowerblue', edgecolor='k')
    axarr[cat_i].bar(ind + 0.2, b_adj_fracs, 0.4, 
                     label=cat[1], color='tomato', edgecolor='k')

    _,p,_,_ = ss.chi2_contingency([[a_fore, b_fore], [a_back, b_back]])

    axarr[cat_i].legend(frameon=False)
    leg = axarr[cat_i].legend()
    leg.set_title("MGP strains", prop={"size": 16})
    if cat_i == 0:
        axarr[cat_i].set_ylabel('Fraction of singletons', fontsize=16)
    
    axarr[cat_i].set_xticks(np.arange(6))
    axarr[cat_i].set_xticklabels([m.replace('>', r'$\to$') for m in mut2idx], fontsize=16)

    sns.despine(ax=axarr[cat_i], top=True, right=True)
        
f.savefig(args.out, bbox_inches='tight')
