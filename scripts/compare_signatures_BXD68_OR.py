import pandas as pd
import scipy.stats as ss
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
import numpy as np
import argparse
from utils import revcomp

plt.rc('font', size=16)

p = argparse.ArgumentParser()
p.add_argument("--annotated_singletons", required=True, 
                    help="""annotated singleton variants in extended BED format.""")
p.add_argument("--out", required=True,
                    help="""name of output plots""")
args = p.parse_args()

# read in BXD singleton data
singleton = pd.read_csv(args.annotated_singletons)

singleton_no_bxd68 = singleton.query('bxd_strain != "4512-JFI-0462_BXD68_RwwJ_phased_possorted_bam"')
singleton_bxd68 = singleton.query('bxd_strain == "4512-JFI-0462_BXD68_RwwJ_phased_possorted_bam"')

group_cols = ['kmer']

# convert to wide-form dataframe
singleton_bxd68_tidy = singleton_bxd68.groupby(group_cols).count().add_suffix('_count').reset_index()

# subset tidy dataframe to relevant columns
group_cols.append('chrom_count')
singleton_bxd68_tidy = singleton_bxd68_tidy[group_cols]

bxd68_total = np.sum(singleton_bxd68_tidy['chrom_count'])

singleton_bxd68_tidy['frac'] = singleton_bxd68_tidy['chrom_count'] / bxd68_total

group_cols = ['kmer', "haplotype_at_qtl"]

# convert to wide-form dataframe
singleton_no_bxd68_tidy = singleton_no_bxd68.groupby(group_cols).count().add_suffix('_count').reset_index()

# subset tidy dataframe to relevant columns
group_cols.append('chrom_count')
singleton_no_bxd68_tidy = singleton_no_bxd68_tidy[group_cols]

# generate subsets of variants in each of two categories, defined
# by the two unique values that the `subset_key` column can take on
subset_0 = singleton_no_bxd68_tidy[singleton_no_bxd68_tidy["haplotype_at_qtl"] == "B"]['chrom_count'].values
subset_1 = singleton_no_bxd68_tidy[singleton_no_bxd68_tidy["haplotype_at_qtl"] == "D"]['chrom_count'].values

# convert counts in the subsets to fractions
subset_0_fracs = subset_0 / np.sum(subset_0)
subset_1_fracs = subset_1 / np.sum(subset_1)

# calculate significant odds-ratio differences between
# the two subsets
pvals = np.ones(subset_0.shape[0], dtype=np.float64)

for i in np.arange(subset_0.shape[0]):
    a_fore, b_fore = subset_0[i], subset_1[i]
    a_back = np.sum(subset_0) - a_fore
    b_back = np.sum(subset_1) - b_fore

    _, p, _, _ = ss.chi2_contingency([[a_fore, b_fore], [a_back, b_back]])

    pvals[i] = p

# convert fractions to log-2 ratios
ratios = subset_1_fracs / subset_0_fracs
log_ratios = np.log2(ratios)

# get a mapping of each mutation type to a corresponding index
uniq_kmers = list(pd.unique(singleton_no_bxd68_tidy['kmer']))
mut2idx = dict(zip(uniq_kmers, range(len(uniq_kmers))))

# get a list of the 6 "base" mutation types in the signature
base_muts = ['>'.join([m.split('>')[0][1], m.split('>')[1][1]]) for m in uniq_kmers]
base_muts = set(list(base_muts))

f, ax = plt.subplots(figsize=(8,8))

sns.set_style('ticks')

ind = np.arange(len(mut2idx))

# get number of significantly enriched mutation types in the BXD
n_sig = np.where(pvals < 0.05 / 96)[0].shape[0]

colors = sns.color_palette('colorblind', n_sig)

sig_counted = 0

for mut in mut2idx:
    bxd68_sub = singleton_bxd68_tidy[singleton_bxd68_tidy['kmer'] == mut]
    # we can only compare spectra for mutations that are in both
    # BXD68 and the rest of the dataset
    if bxd68_sub.shape[0] == 0: continue
    x = bxd68_sub['frac'].values[0]

    idx = mut2idx[mut]

    # manual adjustments so that text annotations look OK
    mut2format = {"TCT>TAT": (10, 30),
                  "TCA>TAA": (-40, 40),
                  "TCC>TAC": (-40, 20),
                  "GCA>GAA": (-42, 30),
                  "GCT>GAT": (5, 50),
                  "CCA>CAA": (-100, -60),
                  "CCT>CAT": (-125, 35)}

    y = log_ratios[mut2idx[mut]]

    edgecolor = "w"
    label = None
    color = "grey"
    if pvals[idx] < 0.05 / 96: 
        edgecolor = 'k'
        label = mut
        color = "firebrick"
        sig_counted += 1

    ax.scatter(x, y, color=color, s=200, edgecolor=edgecolor)

    if pvals[idx] < 0.05 / 96:
        text = mut
        try:
            xytext = mut2format[text]
        except KeyError: continue
        text = text.replace('>', r"$\to$")
        ax.annotate(text,
                    (x, y),
                    xytext=xytext, 
                    arrowprops=dict(facecolor='k',
                    headwidth=0.1, headlength=0.2, lw=0.5),
                    textcoords='offset points', zorder=0)

ax.set_xlabel('Fraction of BXD68 singletons', fontsize=18)
ax.set_ylabel(r'$log_{2}$' + ' ratio of singleton fractions\nin strains with D vs. B haplotypes at QTL', fontsize=18)

sns.despine(ax=ax, top=True, right=True)

f.savefig(args.out, bbox_inches='tight')

