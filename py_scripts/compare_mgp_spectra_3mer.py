import argparse
import scipy.stats as ss
import itertools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from figure_gen_utils import revcomp
from mutation_comparison import mutation_comparison

plt.rc('font', size=18)

p = argparse.ArgumentParser()
p.add_argument("--dumont_xls", required=True,
                    help="""strain-private substitution data from Dumont et al. (2019)""")
p.add_argument("--annotated_singletons", required=True, 
                        help="""annotated singleton variants in extended BED format.""")
p.add_argument("--out", required=True,
                    help="""name of output plot""")
args = p.parse_args()

def get_kmer(row):
    anc, der = row['ancestral'], row['derived']

    kmer_anc = row['5prime_flank'] + anc + row['3prime_flank']
    kmer_der = row['5prime_flank'] + der + row['3prime_flank']

    if anc not in ("C", "A"):
        kmer_anc = revcomp(kmer_anc)
        kmer_der = revcomp(kmer_der)
    return "{}>{}".format(kmer_anc, kmer_der)


# read in Dumont singleton data
dumont = pd.read_excel(args.dumont_xls, 
                       sheet_name="TableS4",
                       header=2)

dumont['kmer'] = dumont.apply(lambda row: get_kmer(row), axis=1)

mut2idx = dict(zip(pd.unique(dumont['kmer']), range(len(pd.unique(dumont['kmer'])))))
idx2mut = {v:k for k,v in mut2idx.items()}

dumont['kmer_idx'] = dumont['kmer'].apply(lambda k: mut2idx[k])

dumont = dumont.sort_values('kmer_idx')

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
cat2strain = {'5/5': withvar_all, '3/5': withvar_some, '2/5': loner, '0/5': novar}

# subset dumont strains by how DBA/2J or C57BL/6J-like they are
a = dumont[dumont['FocalStrain'].isin(withvar_all)]
b = dumont[dumont['FocalStrain'].isin(novar)]

# count total numbers of each kmer mutation type in each group
a_tidy = a.groupby('kmer_idx').count().add_suffix('_count').reset_index()
b_tidy = b.groupby('kmer_idx').count().add_suffix('_count').reset_index()

subset_0 = a_tidy['chr_count'].values
subset_1 = b_tidy['chr_count'].values

# convert sums to log-ratios of fractions of each mutation type
subset_0_fracs = subset_0 / np.sum(subset_0)
subset_1_fracs = subset_1 / np.sum(subset_1)

dumont_log_ratios = np.log2(subset_0_fracs / subset_1_fracs)

# read in the BXD singleton data
singleton = pd.read_csv(args.annotated_singletons)

group_cols = ['kmer', "haplotype_at_qtl"]

# convert to wide-form dataframe
singleton_tidy = singleton.groupby(group_cols).count().add_suffix('_count').reset_index()

# subset tidy dataframe to relevant columns
group_cols.append('chrom_count')
singleton_tidy = singleton_tidy[group_cols]

singleton_tidy['kmer_idx'] = singleton_tidy['kmer'].apply(lambda k: mut2idx[k])
singleton_tidy = singleton_tidy.sort_values('kmer_idx')

# generate subsets of variants in each of two categories, defined
# by the two unique values that the `subset_key` column can take on
subset_0 = singleton_tidy[singleton_tidy["haplotype_at_qtl"] == "B"]['chrom_count'].values
subset_1 = singleton_tidy[singleton_tidy["haplotype_at_qtl"] == "D"]['chrom_count'].values

# make sure both datasets have the same mutation type at each index
muts_0 = singleton_tidy[singleton_tidy["haplotype_at_qtl"] == "B"]['kmer'].values
muts_1 = singleton_tidy[singleton_tidy["haplotype_at_qtl"] == "D"]['kmer'].values

# convert counts in the subsets to log-ratio fractions
subset_0_fracs = subset_0 / np.sum(subset_0)
subset_1_fracs = subset_1 / np.sum(subset_1)

ratios = subset_1_fracs / subset_0_fracs
singleton_log_ratios = np.log2(ratios)

# calculate significant odds-ratio differences between
# the two subsets
pvals = np.ones(subset_0.shape[0], dtype=np.float64)

for i in np.arange(subset_0.shape[0]):
    a_fore, b_fore = subset_0[i], subset_1[i]
    a_back = np.sum(subset_0) - a_fore
    b_back = np.sum(subset_1) - b_fore

    _, p, _, _ = ss.chi2_contingency([[a_fore, b_fore], [a_back, b_back]])

    pvals[i] = p

# make figure
f, ax = plt.subplots(figsize=(10,8))

colors = sns.color_palette('colorblind', 2)

for mut in mut2idx:
    mut_i = mut2idx[mut]
    
    edgecolor = "w"
    label = None
    color = "grey"
    if pvals[mut_i] < 0.05 / 96: 
        edgecolor = 'k'
        label = mut
        color = "firebrick"

    x = singleton_log_ratios[mut_i]
    y = dumont_log_ratios[mut_i]

    ax.scatter(x, 
            y, color=color,
            edgecolor=edgecolor, label=label, s=100)

    # manual adjustments so that text annotations look OK
    mut2format = {"TCT>TAT": (10, 30),
                "TCA>TAA": (-50, -60),
                "TCC>TAC": (-50, 20),
                "GCA>GAA": (-42, 40),
                "GCT>GAT": (40, 30),
                "CCA>CAA": (-40, 60),
                "CCT>CAT": (40, -50)}

    if pvals[mut_i] < 0.05 / 96:
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

ax.set_ylabel('Log-2 ratio of mutation fractions in DBA-like\nvs.C57-like MGP strains')
ax.set_xlabel('Log-2 ratio of mutation fractions in BXDs\nwith D vs. B haplotypes at QTL')
sns.set_style('ticks')
sns.despine(top=True, right=True)
f.savefig(args.out, dpi=300, bbox_inches='tight')