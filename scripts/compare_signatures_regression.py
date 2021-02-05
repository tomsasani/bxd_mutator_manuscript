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
p.add_argument("--out", required=True,
                    help="""name of output plots""")
p.add_argument("--sig_name", required=True,
                    help="""name of COSMIC signature [SBS36_mm10, SBS18_mm10]""")
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
cosmic_components = cosmic[args.sig_name].values

# make figure object

f, ax = plt.subplots(figsize=(6,8))

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

reordered_p = pvals[reordered_idxs][sorted_cosmic_idxs]

edgecolors = ['k' if reordered_p[i] < 0.05 / 96 else 'w' for i in np.arange(reordered_p.shape[0])]

for idx in sorted_cosmic_idxs:

    x = reordered_ratios[sorted_cosmic_idxs][idx]
    y = cosmic_components[sorted_cosmic_idxs][idx]
    c = bar_colors[sorted_cosmic_idxs][idx]
    ec = 'w'
    s = 100
    if reordered_p[idx] < 0.05 / 96: ec = 'k'

    mut2format = {"TCT>TAT": (-40, 20),
                  "TCA>TAA": (-40, 20),
                  "TCC>TAC": (-40, 20),
                  "GCA>GAA": (-42, 50),
                  "GCT>GAT": (5, 20),
                  "CCA>CAA": (-40, -60),
                  "CCT>CAT": (5, -35)}

    if reordered_p[idx] < 0.05 / 96:
        ec = 'k'
        text = pd.unique(cosmic['kmer'])[sorted_cosmic_idxs][idx]
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
    ax.scatter(x, y, edgecolor=ec, s=s, c=c)


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
legend_elements = [Patch(facecolor=c, edgecolor='w', label=c2mut[c]) for c in c2mut]
ax.legend(handles=legend_elements, frameon=False, 
            fontsize=16)

ax.set_ylabel('Fraction of COSMIC {} signature'.format(args.sig_name.split('_')[0]), fontsize=18)
ax.set_xlabel(r'$log_{2}$' + ' ratio of singleton fractions\nin strains with D vs. B haplotypes at QTL', fontsize=18)

sns.despine(ax=ax, top=True, right=True)

f.savefig(args.out, bbox_inches='tight')

