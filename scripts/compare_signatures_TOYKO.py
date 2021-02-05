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
p.add_argument("--ohno_data", required=True,
                    help="""supplementary file from ohno et al. (2014).""")
p.add_argument("--out", required=True,
                    help="""name of output plots""")
args = p.parse_args()

def revcomp(nuc):
    d = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

    return ''.join([d[n] for n in list(nuc)])[::-1] 

def convert_toyko_mutation(sequence):
    mutation = sequence.split('[')[-1].split(']')[0]
    left_flank_1bp = sequence.split('/')[0].split('[')[0][-1]
    right_flank_1bp = sequence.split('/')[-1].split(']')[-1][0]

    #print (mutation, left_flank_1bp, right_flank_1bp)

    anc, der = mutation.split('/')

    kmer_anc = left_flank_1bp + anc + right_flank_1bp 
    kmer_der = left_flank_1bp + der + right_flank_1bp

    if mutation not in ["C/A", "G/T"]: return 'not_CA'

    # reverse complement if necessary
    rc = False
    if mutation[0] == "G":
        rc = True
 
    
    if rc: return "{}>{}".format(revcomp(kmer_anc), revcomp(kmer_der))
    else: return "{}>{}".format(kmer_anc, kmer_der)



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

# read in the TOY-KO signature
toyko = pd.read_excel(args.ohno_data, header=4)

colname = "SEQUENCE  (50+[W/M]+50)"

toyko = toyko[[colname]]
toyko = toyko.fillna(0)
toyko = toyko[toyko[colname] != 0]

toyko['kmer'] = toyko[colname].apply(lambda s: convert_toyko_mutation(s))

toyko = toyko.query('kmer != "not_CA"')

toyko_wide = toyko.groupby('kmer').count().add_suffix('_count').reset_index()

toyko_total = np.sum(toyko_wide[colname + '_count'])

toyko_wide['frac'] = toyko_wide[colname + '_count'] / toyko_total

f, ax = plt.subplots(figsize=(8,8))

sns.set_style('ticks')

ind = np.arange(len(mut2idx))

n_sig = np.where(pvals < 0.05 / 96)[0].shape[0]

colors = sns.color_palette('colorblind', n_sig)

sig_counted = 0

for mut in mut2idx:
    toyko_sub = toyko_wide[toyko_wide['kmer'] == mut]
    if toyko_sub.shape[0] == 0: continue
    toyko_frac = toyko_sub['frac'].values[0]

    bxd_or = log_ratios[mut2idx[mut]]

    bxd_pval = pvals[mut2idx[mut]]

    edgecolor = "w"
    label = None
    color = "grey"
    if bxd_pval < 0.05 / 96: 
        edgecolor = 'k'
        label = mut
        color = colors[sig_counted]
        sig_counted += 1

    ax.scatter(bxd_or, toyko_frac, c=color, s=150, edgecolor=edgecolor, label=label)



ax.legend()

# create custom legend
#legend_elements = [Patch(facecolor=c, edgecolor='w', label=c2mut[c]) for c in c2mut]
#ax.legend(handles=legend_elements, frameon=False, 
#            fontsize=16)

ax.set_ylabel('Fraction of de novo germline\nmutations in TOY-KO mice', fontsize=18)
ax.set_xlabel(r'$log_{2}$' + ' ratio of singleton fractions\nin strains with D vs. B haplotypes at QTL', fontsize=18)

sns.despine(ax=ax, top=True, right=True)

f.savefig(args.out, bbox_inches='tight')

