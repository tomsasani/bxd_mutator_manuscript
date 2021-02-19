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
                    help="""name of output plot""")
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

singleton = singleton.query('bxd_strain == "4512-JFI-0462_BXD68_RwwJ_phased_possorted_bam"')

group_cols = ['kmer']

# convert to wide-form dataframe
singleton_tidy = singleton.groupby(group_cols).count().add_suffix('_count').reset_index()

# subset tidy dataframe to relevant columns
group_cols.append('chrom_count')
singleton_tidy = singleton_tidy[group_cols]

bxd68_total = np.sum(singleton_tidy['chrom_count'])

singleton_tidy['frac'] = singleton_tidy['chrom_count'] / bxd68_total

# get a mapping of each mutation type to a corresponding index
uniq_kmers = list(pd.unique(singleton_tidy['kmer']))
mut2idx = dict(zip(uniq_kmers, range(len(uniq_kmers))))

# get a list of the 6 "base" mutation types in the signature
base_muts = ['>'.join([m.split('>')[0][1], m.split('>')[1][1]]) for m in uniq_kmers]
base_muts = list(set(base_muts))

# read in the COSMIC signature
cosmic = pd.read_csv(args.cosmic_signature)

# convert COSMIC mutation notation to match mine
cosmic['kmer'] = cosmic.apply(convert_cosmic_mutation, axis=1)

# sort COSMIC mutations with respect to the order of mutations in
# the `log_ratios` numpy array
cosmic['kmer_idx'] = cosmic['kmer'].apply(lambda k: mut2idx[k] if k in mut2idx else 'NA')

# remove kmers for which we don't have data in BXD68
cosmic = cosmic[cosmic['kmer_idx'] != "NA"]

# make figure object
f, ax = plt.subplots(figsize=(6,8))

sns.set_style('ticks')

colors = sns.color_palette('colorblind', len(base_muts))
mut2c = dict(zip(base_muts, colors))

for mut in mut2idx:

    idx = mut2idx[mut]

    bxd68_sub = singleton_tidy[singleton_tidy['kmer'] == mut]
    x = bxd68_sub['frac'].values[0]    
    
    y = cosmic[cosmic['kmer'] == mut][args.sig_name].values[0]

    nuc_a = mut.split('>')[0][1]
    nuc_b = mut.split('>')[1][1]
    base_mut = "{}>{}".format(nuc_a, nuc_b)

    c = mut2c[base_mut]

    s = 100

    # manual adjustments so that text annotations look OK
    mut2format = {"TCT>TAT": (-40, 20),
                  "TCA>TAA": (-40, 20),
                  "TCC>TAC": (-40, 20),
                  "GCA>GAA": (-42, 50),
                  "GCT>GAT": (5, 20),
                  "CCA>CAA": (-40, -60),
                  "CCT>CAT": (5, -35)}

    ax.scatter(x, y, edgecolor='w', s=s, color=c)

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

# create custom legend
legend_elements = [Patch(facecolor=mut2c[mut], edgecolor='w', label=mut) for mut in mut2c]
ax.legend(handles=legend_elements, frameon=False, 
            fontsize=16)

ax.set_ylabel('Fraction of COSMIC {} signature'.format(args.sig_name.split('_')[0]), fontsize=18)
ax.set_xlabel(r'Fraction of BXD68 singletons', fontsize=18)

sns.despine(ax=ax, top=True, right=True)

f.savefig(args.out, bbox_inches='tight')

