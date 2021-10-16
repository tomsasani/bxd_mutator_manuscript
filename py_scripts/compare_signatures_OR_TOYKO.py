import pandas as pd
import scipy.stats as ss
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
import numpy as np
import argparse
from figure_gen_utils import revcomp

plt.rc('font', size=16)

p = argparse.ArgumentParser()
p.add_argument("--annotated_singletons", required=True, 
                    help="""annotated singleton variants in extended BED format.""")
p.add_argument("--ohno_data", required=True,
                    help="""supplementary file from ohno et al. (2014).""")
p.add_argument("--out", required=True,
                    help="""name of output plots""")
args = p.parse_args()

def convert_toyko_mutation(sequence):
    """
    convert TOY-KO mutations (reported as a string of 50
    upstream nucleotides plus the mutation plus a string of 50
    downstream nucleotides) to notation that matches the BXD data
    """
    mutation = sequence.split('[')[-1].split(']')[0]
    left_flank_1bp = sequence.split('/')[0].split('[')[0][-1]
    right_flank_1bp = sequence.split('/')[-1].split(']')[-1][0]

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

# read in BXD singleton data
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
muts_b = singleton_tidy[singleton_tidy["haplotype_at_qtl"] == "B"]["kmer"].values
muts_d = singleton_tidy[singleton_tidy["haplotype_at_qtl"] == "D"]["kmer"].values

assert np.array_equal(muts_b, muts_d)

# get a mapping of each mutation type to a corresponding index
mut2idx = dict(zip(muts_b, range(len(muts_b))))

# convert counts in the subsets to fractions
subset_0_fracs = subset_0 / np.sum(subset_0)
subset_1_fracs = subset_1 / np.sum(subset_1)

# calculate significant odds-ratio differences between
# the two subsets
pvals = np.ones(subset_0.shape[0], dtype=np.float64)

for mut in mut2idx:
    i = mut2idx[mut]
    a_fore, b_fore = subset_0[i], subset_1[i]
    a_back = np.sum(subset_0) - a_fore
    b_back = np.sum(subset_1) - b_fore

    _, p, _, _ = ss.chi2_contingency([[a_fore, b_fore], [a_back, b_back]])

    pvals[i] = p

# convert fractions to log-2 ratios
ratios = subset_1_fracs / subset_0_fracs
log_ratios = np.log2(ratios)

# get a list of the 6 "base" mutation types in the signature
base_muts = ['>'.join([m.split('>')[0][1], m.split('>')[1][1]]) for m in muts_b]
base_muts = base_muts[::16]

# read in the TOY-KO signature
toyko = pd.read_excel(args.ohno_data, header=4)

# filter the TOY-KO data to only include the 
# individual mutations observed, plus or minus 50bp of context
colname = "SEQUENCE  (50+[W/M]+50)"
toyko = toyko[[colname]]
toyko = toyko.fillna(0)
toyko = toyko[toyko[colname] != 0]

# convert the raw sequence to a 3-mer mutation
toyko['kmer'] = toyko[colname].apply(lambda s: convert_toyko_mutation(s))

# only consider the C>A mutation types (99% of the data)
toyko = toyko.query('kmer != "not_CA"')

# get the relative frequencies of each 3-mer C>A mutation
toyko_wide = toyko.groupby('kmer').count().add_suffix('_count').reset_index()
toyko_total = np.sum(toyko_wide[colname + '_count'])
toyko_wide['frac'] = toyko_wide[colname + '_count'] / toyko_total

f, ax = plt.subplots(figsize=(8,8))

sns.set_style('ticks')

ca_mut2idx = {k:v for k,v in mut2idx.items() if k[1] == "C" and k[5] == "A"}
ind = np.arange(len(ca_mut2idx))

# get number of significantly enriched mutation types in the BXD
n_sig = np.where(pvals < 0.05 / 96)[0].shape[0]

sig_counted = 0

a, b = [], []

for mut in ca_mut2idx:
    toyko_sub = toyko_wide[toyko_wide['kmer'] == mut]
    if toyko_sub.shape[0] == 0: 
        toyko_frac = 0.
    else: toyko_frac = toyko_sub['frac'].values[0]

    idx = mut2idx[mut]

    # manual adjustments so that text annotations look OK
    mut2format = {"TCT>TAT": (10, 30),
                  "TCA>TAA": (-40, 40),
                  "TCC>TAC": (-40, 20),
                  "GCA>GAA": (-42, -50),
                  "GCT>GAT": (-60, 50),
                  "CCA>CAA": (-40, -60),
                  "CCT>CAT": (25, -35)}

    x = log_ratios[mut2idx[mut]]
    y = toyko_frac

    a.append(x)
    b.append(y)

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

ax.set_ylabel('Fraction of de novo germline\nmutations in TOY-KO mice', fontsize=18)
ax.set_xlabel('Log-2 ratio of ' + r'C$\to$A' + ' singleton fractions\nin strains with D vs. B haplotypes at QTL', fontsize=18)

sns.despine(ax=ax, top=True, right=True)
print (ss.spearmanr(a, b))
slope, intercept, r_value, p_value, std_err = ss.linregress(a, b)
print (p_value)
f.savefig(args.out, bbox_inches='tight')

