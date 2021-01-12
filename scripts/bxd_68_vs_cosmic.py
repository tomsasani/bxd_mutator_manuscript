import pandas as pd
import scipy.stats as ss
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
import numpy as np
import argparse

plt.rc('font', size=16)

p = argparse.ArgumentParser()
p.add_argument("--annotated_singletons", required=True, 
                    help="""annotated singleton variants in extended BED format.""")
p.add_argument("--cosmic_signature", required=True,
                    help="""file containing mutation spectrum corresponding to SBS36.""")
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

singleton = singleton.query('haplotype_at_qtl == 1')
#singleton = singleton.query('bxd_strain == "4512-JFI-0462_BXD68_RwwJ_phased_possorted_bam"')

group_cols = ['kmer']

# convert to wide-form dataframe
singleton_tidy = singleton.groupby(group_cols).count().add_suffix('_count').reset_index()

# subset tidy dataframe to relevant columns
group_cols.append('chrom_count')
singleton_tidy = singleton_tidy[group_cols]

# read in the COSMIC signature
cosmic = pd.read_csv(args.cosmic_signature)

# convert COSMIC mutation notation to match mine
cosmic['kmer'] = cosmic.apply(convert_cosmic_mutation, axis=1)

# get a mapping of each mutation type to a corresponding index
uniq_kmers = list(pd.unique(cosmic['kmer']))
mut2idx = dict(zip(uniq_kmers, range(len(uniq_kmers))))

# get a list of the 6 "base" mutation types in the signature
base_muts = ['>'.join([m.split('>')[0][1], m.split('>')[1][1]]) for m in mut2idx]
base_muts = base_muts[::16]

# get the components of the COSMIC mutation signature
cosmic = cosmic.sort_values('Type')
cosmic_components = cosmic['SBS36_mm10'].values

# make figure object

f, ax = plt.subplots(figsize=(6,8))

sns.set_style('ticks')

colors = sns.color_palette('colorblind', 6)
bar_colors = np.repeat(colors, 16, axis=0)

ind = np.arange(len(mut2idx))

# reorder ratios from singletons with respect to the order
# of mutations in the COSMIC signature
reordered_idxs = np.array([mut2idx[m] for m in pd.unique(cosmic['kmer'])])

singleton_complete = np.zeros(ind.shape[0])
for m in pd.unique(cosmic['kmer']):
    idx = mut2idx[m]
    s_count = singleton_tidy[singleton_tidy['kmer'] == m]['chrom_count'].values
    if s_count.shape[0] == 0: continue
    singleton_complete[idx] = s_count
singleton_fracs = singleton_complete / np.sum(singleton_complete)

reordered_fracs = singleton_fracs[reordered_idxs]


# sort the COSMIC signature in decreasing order of each kmer's
# contribution to the signature
sorted_cosmic_idxs = np.argsort(cosmic_components)

singleton_ranks = np.argsort(reordered_fracs)
cosmic_ranks = np.argsort(cosmic_components)

print (ss.kendalltau(singleton_ranks, cosmic_ranks))

ax.scatter(reordered_fracs[sorted_cosmic_idxs], 
           cosmic_components[sorted_cosmic_idxs], 
           c=bar_colors[sorted_cosmic_idxs], edgecolor='w', s=100)

# map colors to the base mutation type they correspond to
base_muts = []
for m in pd.unique(cosmic['kmer']):
    nuc_a = m.split('>')[0][1]
    nuc_b = m.split('>')[1][1]

    base_mut = r'$\to$'.join([nuc_a, nuc_b])

    if base_mut in base_muts: continue
    base_muts.append(base_mut)
c2mut = dict(zip(colors, base_muts))

# create custom legend
legend_elements = [Patch(facecolor=c, edgecolor='w', label=c2mut[c]) for c in c2mut]
ax.legend(handles=legend_elements, frameon=False, 
            fontsize=16)

ax.set_ylabel('Fraction of COSMIC SBS36 signature', fontsize=18)
ax.set_xlabel(r'$log_{2}$' + ' ratio of singleton fractions\nin strains with D2 vs. B6 haplotypes at QTL', fontsize=18)
f.savefig('oo.eps', bbox_inches='tight')

