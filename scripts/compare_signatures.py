import pandas as pd
import scipy.stats as ss
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
import numpy as np
import argparse

plt.rc('font', size=12)

p = argparse.ArgumentParser()
p.add_argument("--annotated_singletons", required=True, 
                    help="""annotated singleton variants in extended BED format.""")
p.add_argument("--cosmic_signature", required=True,
                    help="""file containing mutation spectrum corresponding to SBS36.""")
p.add_argument("--out", required=True,
                    help="""name of output plots""")
p.add_argument("-orientation", required=False, default="h",
                    help="""make plot horizontal or vertical [h,v]""")
args = p.parse_args()

def revcomp(nuc):
    d = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

    return ''.join([d[n] for n in list(nuc)])[::-1] 

def convert_cosmic_mutation(row):
    context = row['Subtype']
    mutation = row['Type']

    changed_from = mutation.split('>')[0]
    changed_to = mutation.split('>')[1]

    e_kmer = context[0] + changed_to + context[-1]

    if changed_from == "T":
        context = revcomp(context)
        e_kmer = revcomp(e_kmer)

    return context + '>' + e_kmer


singleton = pd.read_csv(args.annotated_singletons)

group_cols = ['kmer', "haplotype_at_qtl"]

# convert to wide-form dataframe
singleton_tidy = singleton.groupby(group_cols).count().add_suffix('_count').reset_index()

# subset tidy dataframe to relevant columns
group_cols.append('chrom_count')
singleton_tidy = singleton_tidy[group_cols]

# generate subsets of variants in each of two categories, defined
# by the two unique values that the `subset_key` column can take on
subset_0 = singleton_tidy[singleton_tidy["haplotype_at_qtl"] == 0]['chrom_count'].values
subset_1 = singleton_tidy[singleton_tidy["haplotype_at_qtl"] == 1]['chrom_count'].values

# convert counts in the subsets to fractions
subset_0_fracs = subset_0 / np.sum(subset_0)
subset_1_fracs = subset_1 / np.sum(subset_1)

# convert fractions to log-2 ratios
ratios = subset_1_fracs / subset_0_fracs
log_ratios = np.log2(ratios)

# get a mapping of each mutation type to a corresponding index
uniq_kmers = list(pd.unique(singleton_tidy['kmer']))
mut2idx = dict(zip(uniq_kmers, range(len(uniq_kmers))))

# get a list of the 6 "base" mutation types in the signature
base_muts = ['>'.join([m.split('>')[0][1], m.split('>')[1][1]]) for m in uniq_kmers]
base_muts = base_muts[::16]

# read in the COSMIC signature
cosmic = pd.read_csv(args.cosmic_signature)

# convert COSMIC mutation notation to match mine
cosmic['kmer'] = cosmic.apply(convert_cosmic_mutation, axis=1)

# get the components of the COSMIC mutation signature
cosmic = cosmic.sort_values('Type')
cosmic_components = cosmic['SBS36_mm10'].values

# make figure object

if args.orientation == "v":
    f, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 10))
elif args.orientation == "h":
    f, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

sns.set_style('ticks')

colors = sns.color_palette('colorblind', 6)
bar_colors = np.repeat(colors, 16, axis=0)

ind = np.arange(len(mut2idx))

# reorder ratios from singletons with respect to the order
# of mutations in the COSMIC signature
reordered_idxs = np.array([mut2idx[m] for m in pd.unique(cosmic['kmer'])])
reordered_ratios = log_ratios[reordered_idxs]



# sort the COSMIC signature in decreasing order of each kmer's
# contribution to the signature
sorted_cosmic_idxs = np.argsort(cosmic_components)

# plot the COSMIC components and my singleton enrichments in
# the same order

if args.orientation == "v":
    ax1.barh(ind, cosmic_components[sorted_cosmic_idxs], 1, 
                color=bar_colors[sorted_cosmic_idxs], edgecolor='k')
    ax2.barh(ind, reordered_ratios[sorted_cosmic_idxs], 1, 
                color=bar_colors[sorted_cosmic_idxs], edgecolor='k')

elif args.orientation == "h":
    ax1.bar(ind, cosmic_components[sorted_cosmic_idxs], 1, 
                color=bar_colors[sorted_cosmic_idxs], edgecolor='k')
    ax2.bar(ind, reordered_ratios[sorted_cosmic_idxs], 1, 
                color=bar_colors[sorted_cosmic_idxs], edgecolor='k')

# map colors to the base mutation type they correspond to
base_muts = []
for m in pd.unique(cosmic['kmer'])[sorted_cosmic_idxs][::-1]:
    nuc_a = m.split('>')[0][1]
    nuc_b = m.split('>')[1][1]

    base_mut = r'$\to$'.join([nuc_a, nuc_b])

    if base_mut in base_muts: continue
    base_muts.append(base_mut)
c2mut = dict(zip(colors, base_muts))

# create custom legend
legend_elements = [Patch(facecolor=c, edgecolor='k', label=c2mut[c]) for c in c2mut]
ax1.legend(handles=legend_elements, frameon=False, 
            fontsize=14, loc="center")

# add extra plot stuff
if args.orientation == "v":
    ax1.set_yticks(ind)
    ax1.set_yticklabels([m.replace('>', r'$\to$') for m in pd.unique(cosmic['kmer'])[sorted_cosmic_idxs]], fontsize=6)
    ax1.set_xlabel('Fraction of COSMIC\nSBS36 signature', fontsize=14)
    ax2.set_xlabel(r'$log_{2}$' + ' ratio of singleton fractions\nin strains with D2 vs. B6 haplotypes at QTL', fontsize=14)
    ax2.set_yticks([])
elif args.orientation == "h":
    ax2.set_xticks(ind)
    ax2.set_xticklabels([m.replace('>', r'$\to$') for m in pd.unique(cosmic['kmer'])[sorted_cosmic_idxs]], fontsize=6, rotation=90)
    ax1.set_ylabel('Fraction of COSMIC\nSBS36 signature', fontsize=14)
    ax2.set_ylabel(r'$log_{2}$' + ' ratio of singleton fractions\nin strains with D2 vs. B6\nhaplotypes at QTL', fontsize=14)
    ax1.set_xticks([])

sns.despine(ax=ax1, top=True, right=True)
sns.despine(ax=ax2, top=True, right=True)

f.savefig(args.out, bbox_inches='tight')

