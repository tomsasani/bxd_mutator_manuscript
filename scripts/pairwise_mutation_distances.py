import pandas as pd
import itertools
import matplotlib.pyplot as plt
from collections import defaultdict
import seaborn as sns
import numpy as np
import argparse

plt.rc('font', size=12)

p = argparse.ArgumentParser()
p.add_argument("--annotated_singletons", required=True, 
                    help="""annotated singleton variants in extended BED format.""")
p.add_argument("--out", required=True,
                    help="""name of output plot""")
args = p.parse_args()

# read in singletons
singleton = pd.read_csv(args.annotated_singletons)

# map strains to corresponding indices
strains = pd.unique(singleton['bxd_strain_conv'])
strain2idx = dict(zip(strains, range(len(strains))))

chroms = pd.unique(singleton['chrom'])

# map mutations to corresponding indices
mut2idx = dict(zip(["C>A", "C>T", "C>G", 
                    "A>T", "A>C", "A>G", "CpG>TpG"], range(7)))

# store the total numbers of all singletons and of all
# dicnucleotide mutations per strain
total_muts_per_strain = np.zeros(len(strains))
total_mnms_per_strain = np.zeros(len(strains))

# loop over every strain
for strain in strains:
    strain_idx = strain2idx[strain]
    # within each strain, loop over every chromosome
    for chrom in chroms:
        # subset that strain's singletons to just the ones that
        # occur on this chromosome
        df_sub = singleton[ (singleton['bxd_strain_conv'] == strain) & \
                            (singleton['chrom'] == chrom) ]

        # get an array of the start positions of every singleton in
        # this strain on this chromosome
        starts = df_sub['start'].values

        # increment this strain's total count of MNMs by the length
        # of that array
        total_muts_per_strain[strain_idx] += starts.shape[0]

        # get the distances between each variant and the variant before it
        subsequent_vdist = starts[1:] - starts[:-1]

        # we're only interested in the sites where there are dinucleotide 
        # mutations, and therefore only in the places where the distance
        # between subsequent mutations is exactly 1
        mnm_locations = np.where(subsequent_vdist == 1)[0]

        # if there aren't any MNMs, move on
        if mnm_locations.shape[0] == 0: continue

        # otherwise, increment the strain's count of MNMs
        else:
            total_mnms_per_strain[strain_idx] += mnm_locations.shape[0]

# get the fraction of each strain's singletons that are MNMs
mnm_fracs = total_mnms_per_strain / total_muts_per_strain

# sort strains in order of increasing MNM fractions
sorted_mnm_idxs = np.argsort(mnm_fracs)
sorted_mnm_fracs = mnm_fracs[sorted_mnm_idxs]

ind = np.arange(len(strains))

colors = np.repeat('cornflowerblue', len(strains))

# get the index of the strain we're interested in
bxd68_idx = np.where(strains[sorted_mnm_idxs] == "BXD68")[0]

colors[bxd68_idx] = 'firebrick'

f, ax = plt.subplots(figsize=(8, 16))

sns.set_style('ticks')

ax.barh(ind, sorted_mnm_fracs, 1, edgecolor='k', color=colors)

ax.set_yticks(ind)
ax.set_xlabel("Fraction of singletons that are\nmulti-nucleotide mutations")
ax.set_ylabel("BXD RILs")
ax.set_yticklabels(strains[sorted_mnm_idxs], fontsize=6)

sns.despine(ax=ax, top=True, right=True)

f.savefig(args.out, bbox_inches='tight')

