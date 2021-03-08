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

f, ax = plt.subplots(figsize=(6,8))

sns.set_style('ticks')

ind = np.arange(len(mut2idx))

for mut in mut2idx:
    toyko_sub = toyko_wide[toyko_wide['kmer'] == mut]
    if toyko_sub.shape[0] == 0: continue
    toyko_frac = toyko_sub['frac'].values[0]

    bxd68_sub = singleton_tidy[singleton_tidy['kmer'] == mut]
    bxd68_frac = bxd68_sub['frac'].values[0]

    idx = mut2idx[mut]

    # manual adjustments so that text annotations look OK
    mut2format = {"TCT>TAT": (10, 30),
                  "TCA>TAA": (-40, 40),
                  "TCC>TAC": (-80, 30),
                  "GCA>GAA": (-42, -50),
                  "GCT>GAT": (-60, 50),
                  "CCA>CAA": (-40, -60),
                  "CCT>CAT": (25, -35)}

    x = bxd68_frac
    y = toyko_frac

    edgecolor = "w"
    label = None
    color = "grey"

    ax.scatter(x, y, color=color, s=200, edgecolor=edgecolor)

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
ax.set_xlabel('Fraction of BXD68 singletons', fontsize=18)

sns.despine(ax=ax, top=True, right=True)

f.savefig(args.out, bbox_inches='tight')

