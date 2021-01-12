import argparse
import scipy.stats as ss
import itertools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

p = argparse.ArgumentParser()
p.add_argument("--dumont_xls")
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

dba_idx, c57_idx = smp2idx["DBA_2J"], smp2idx["C57BL_6NJ"]

a_count, b_count = counts[dba_idx, :], counts[c57_idx, :]
a_denom, b_denom = denoms[dba_idx, :], denoms[c57_idx, :]

nuc2denom = {'A':0, 'C':1, 'G':2, 'T':3}

nuc2comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

ind = np.arange(len(muts))

a_adj, b_adj = [], []

# for every mutation type, perform a chi-square test of
# enrichment of that mutation type's counts in subset 0 vs.
# subset 1
for i,mut in enumerate(muts):
    base_nuc = mut.split('>')[0]

    # the "foreground" in the Chi-square test is the sum of
    # the counts of X>Y mutations across strains in subset 0,
    # and the same for subset 1 
    a_fore = a_count[i]
    b_fore = b_count[i]

    # the "background" is the sum of callable base pairs that are
    # X nucleotides across strains in subset 0, and same for subset 1.
    # we add to the "background" the sum of callable base pairs that
    # are the complement of X, as well
    a_denom_ = a_denom[nuc2denom[base_nuc]]
    a_denom_ += a_denom[nuc2denom[nuc2comp[base_nuc]]]
    b_denom_ = b_denom[nuc2denom[base_nuc]]
    b_denom_ += b_denom[nuc2denom[nuc2comp[base_nuc]]]

    # adjustments from michael's paper
    if a_denom_ > b_denom_: 
        adj = b_denom_ / a_denom_
        a_fore = a_fore * adj
    elif b_denom_ > a_denom_:
        adj = a_denom_ / b_denom_
        b_fore = b_fore * adj

    a_adj.append(a_fore)
    b_adj.append(b_fore)

for mut in mut2idx:
    mut_i = mut2idx[mut]
    a_fore = a_adj[mut_i]
    a_back = sum(a_adj) - a_fore
    b_fore = b_adj[mut_i]
    b_back = sum(b_adj) - b_fore

    _, p, _, _ = ss.chi2_contingency([[a_fore, b_fore],[a_back, b_back]])
    fisher_odds, p = ss.fisher_exact([[a_fore, b_fore],[a_back, b_back]])

    print (mut, p, fisher_odds, p)

