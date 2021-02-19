import pandas as pd
import scipy.stats as ss
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
import numpy as np
import argparse
from utils import revcomp

def convert_cosmic_mutation(row):
    """
    convert cosmic mutation notation to 3-mer 
    mutation notation so that it matches the BXD data

    >>> convert_cosmic_mutation("C>G,ACA")
    "ACA>AGA"
    >>> convert_cosmic_mutation("T>C,CTG")
    "CAG>CGG"
    """
    
    context = row['Subtype']
    mutation = row['Type']

    changed_from = mutation.split('>')[0]
    changed_to = mutation.split('>')[1]

    e_kmer = context[0] + changed_to + context[-1]

    if changed_from == "T":
        context = revcomp(context)
        e_kmer = revcomp(e_kmer)

    return context + '>' + e_kmer


p = argparse.ArgumentParser()
p.add_argument("--annotated_singletons", required=True, 
                    help="""annotated singleton variants in extended BED format.""")
p.add_argument("--cosmic_signature", required=True,
                    help="""file containing mutation spectrum corresponding to SBS36.""")
p.add_argument("--out", required=True,
                    help="""name of output plots""")
p.add_argument("--sig_name", required=True,
                    help="""name of COSMIC signature [SBS36_mm10, SBS18_mm10]""")
args = p.parse_args()

plt.rc('font', size=16)

singleton = pd.read_csv(args.annotated_singletons)

group_cols = ['kmer', "haplotype_at_qtl"]

# convert to wide-form dataframe
singleton_tidy = singleton.groupby(group_cols).count().add_suffix('_count').reset_index()

# subset tidy dataframe to relevant columns
group_cols.append('chrom_count')
singleton_tidy = singleton_tidy[group_cols]

# generate subsets of variants in each of two categories, defined
# by the two unique values that the `subset_key` column can take on
subset_0 = singleton_tidy[singleton_tidy["haplotype_at_qtl"] == "B"]['chrom_count'].values
subset_1 = singleton_tidy[singleton_tidy["haplotype_at_qtl"] == "D"]['chrom_count'].values

# make sure both datasets have the same mutation type at each index
muts_0 = singleton_tidy[singleton_tidy["haplotype_at_qtl"] == "B"]['kmer'].values
muts_1 = singleton_tidy[singleton_tidy["haplotype_at_qtl"] == "D"]['kmer'].values

assert np.sum(muts_0 == muts_1) == muts_0.shape[0]

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
#uniq_kmers = list(pd.unique(singleton_tidy['kmer']))
mut2idx = dict(zip(muts_0, range(len(muts_0))))

# get a list of the 6 "base" mutation types in the signature
base_muts = list(set(['>'.join([m.split('>')[0][1], m.split('>')[1][1]]) for m in muts_0]))

# read in the COSMIC signature
cosmic = pd.read_csv(args.cosmic_signature)

# convert COSMIC mutation notation to match mine
cosmic['kmer'] = cosmic.apply(convert_cosmic_mutation, axis=1)

# sort COSMIC mutations with respect to the order of mutations in
# the `log_ratios` numpy array
cosmic['kmer_idx'] = cosmic['kmer'].apply(lambda k: mut2idx[k])
cosmic = cosmic.sort_values('kmer_idx')

# get the components of the COSMIC mutation signature
cosmic_components = cosmic[args.sig_name].values

# make figure object
f, ax = plt.subplots(figsize=(6,8))

sns.set_style('ticks')

colors = sns.color_palette('colorblind', len(base_muts))
mut2c = dict(zip(base_muts, colors))

edgecolors = ['k' if pvals[i] < 0.05 / 96 else 'w' 
                    for i in np.arange(pvals.shape[0])]

for mut in mut2idx:

    idx = mut2idx[mut]

    x = log_ratios[idx]
    y = cosmic_components[idx]

    nuc_a = mut.split('>')[0][1]
    nuc_b = mut.split('>')[1][1]
    base_mut = "{}>{}".format(nuc_a, nuc_b)
    c = mut2c[base_mut]

    ec = 'w'
    s = 100
    if pvals[idx] < 0.05 / 96: ec = 'k'

    # manual adjustments so that text annotations look OK
    mut2format = {"TCT>TAT": (-40, 20),
                  "TCA>TAA": (-40, 20),
                  "TCC>TAC": (-40, 20),
                  "GCA>GAA": (-42, 50),
                  "GCT>GAT": (5, 20),
                  "CCA>CAA": (-40, -60),
                  "CCT>CAT": (5, -35)}

    if pvals[idx] < 0.05 / 96:
        ec = 'k'
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
    ax.scatter(x, y, edgecolor=ec, s=s, color=c)

# create custom legend
legend_elements = [Patch(facecolor=mut2c[mut], edgecolor='w', label=mut) for mut in mut2c]
ax.legend(handles=legend_elements, frameon=False, 
            fontsize=16)

ax.set_ylabel('Fraction of COSMIC {} signature'.format(args.sig_name.split('_')[0]), fontsize=18)
ax.set_xlabel('Log-2 ratio of singleton fractions\nin strains with D vs. B haplotypes at QTL', fontsize=18)

sns.despine(ax=ax, top=True, right=True)

f.savefig(args.out, bbox_inches='tight')

