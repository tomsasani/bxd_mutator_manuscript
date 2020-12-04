import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
import numpy as np
import argparse

plt.rc('font', size=12)

p = argparse.ArgumentParser()
p.add_argument("--cosmic_signature")
p.add_argument("--annotated_singletons")
p.add_argument("--out")
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

tidy_cols = ['kmer', "haplotype_at_qtl"]

# convert to tidy dataframe
singleton_tidy = singleton.groupby(tidy_cols).count().add_suffix('_count').reset_index()

# subset tidy dataframe to relevant columns
tidy_cols.append('chrom_count')

singleton_tidy = singleton_tidy[tidy_cols]

subset_0 = singleton_tidy[singleton_tidy["haplotype_at_qtl"] == 0]['chrom_count'].values
subset_1 = singleton_tidy[singleton_tidy["haplotype_at_qtl"] == 1]['chrom_count'].values

subset_0_fracs = subset_0 / np.sum(subset_0)
subset_1_fracs = subset_1 / np.sum(subset_1)

ratios = subset_1_fracs / subset_0_fracs
log_ratios = np.log2(ratios)

# get a mapping of each mutation type to a corresponding index
uniq_kmers = list(pd.unique(singleton_tidy['kmer']))
mut2idx = dict(zip(uniq_kmers, range(len(uniq_kmers))))

base_muts = ['>'.join([m.split('>')[0][1], m.split('>')[1][1]]) for m in mut2idx]
base_muts = base_muts[::16]

cosmic = pd.read_csv(args.cosmic_signature)

cosmic['kmer'] = cosmic.apply(convert_cosmic_mutation, axis=1)

cosmic = cosmic.sort_values('Type')

cosmic_components = cosmic['SBS36_mm10'].values

f, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 12), sharex=False, sharey=False)

sns.set_style('ticks')

colors = sns.color_palette('colorblind', 6)

bar_colors = np.repeat(colors, 16, axis=0)

ind = np.arange(len(mut2idx))

# reorder ratios from singletons with respect to the order
# of mutations in the COSMIC signature
reordered_idxs = np.array([mut2idx[m] for m in pd.unique(cosmic['kmer'])])

reordered_ratios = log_ratios[reordered_idxs]

sorted_cosmic_idxs = np.argsort(cosmic_components)

ax1.barh(ind, cosmic_components[sorted_cosmic_idxs], 1, 
            color=bar_colors[sorted_cosmic_idxs], edgecolor='k')
ax2.barh(ind, reordered_ratios[sorted_cosmic_idxs], 1, 
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
            fontsize=20, loc="center")

ax1.set_yticks(ind)
ax1.set_yticklabels([m.replace('>', r'$\to$') for m in pd.unique(cosmic['kmer'])[sorted_cosmic_idxs]], fontsize=6)

ax1.set_xlabel('Fraction of COSMIC\nSBS36 signature', fontsize=16)
ax2.set_xlabel(r'$log_{2}$' + ' ratio of singleton fractions\nin strains with D2 vs. B6 haplotypes at QTL', fontsize=16)

ax2.set_yticks([])

sns.despine(ax=ax1, top=True, right=True)
sns.despine(ax=ax2, top=True, right=True)

f.savefig(args.out, bbox_inches='tight')

