import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import scipy.stats as ss
import numpy as np
import seaborn as sns
import pandas as pd
import numpy as np
import math
import argparse
from adj_pvalues import adj_pvalues

p = argparse.ArgumentParser()
p.add_argument("--annotated_singletons", required=True, 
                    help="""annotated singleton variants in extended BED format.""")
p.add_argument("--annotated_dumont", required=True, 
                    help="""annotated singleton variants in extended BED format.""")
p.add_argument("--out", required=True,
                    help="""name of output plot""")
args = p.parse_args()

def to_base_mut(kmer):
    a, b = kmer[1], kmer[5]

    return '{}>{}'.format(a, b)

# read in singletons
singleton = pd.read_csv(args.annotated_singletons)
dumont = pd.read_csv(args.annotated_dumont)

dumont = dumont.replace({"DBA_2J":1, "C57BL_6NJ":0})

singleton['base_mut'] = singleton['kmer'].apply(lambda k: to_base_mut(k))
dumont['base_mut'] = dumont['kmer'].apply(lambda k: to_base_mut(k))
    
# get a mapping of each mutation type to a corresponding index
uniq_kmers = list(pd.unique(singleton['base_mut']))
mut2idx = dict(zip(uniq_kmers, range(len(uniq_kmers))))

x, y = None, None

for data, subset_key, count_key in zip([singleton, dumont], ['haplotype_at_qtl', 'strain'], ['chrom_count', 'value_count']):

    group_cols = ['base_mut', subset_key]

    # convert to wide-form dataframe, grouped by kmer
    df_wide = data.groupby(group_cols).count().add_suffix('_count').reset_index()

    # subset dataframe to relevant columns
    group_cols.append(count_key)
    df_wide = df_wide[group_cols]

    subset_0 = np.zeros(len(mut2idx), dtype=np.int64)
    subset_1 = np.zeros(len(mut2idx), dtype=np.int64)

    odds = np.zeros(len(mut2idx))

    for mut in mut2idx:
        mut_count_a = df_wide[(df_wide[subset_key] == 0) & (df_wide['base_mut'] == mut)][count_key].values
        mut_count_b = df_wide[(df_wide[subset_key] == 1) & (df_wide['base_mut'] == mut)][count_key].values

        a_fore, b_fore = 0, 0
        if mut_count_a.shape[0] > 0: a_fore = mut_count_a[0]
        if mut_count_b.shape[0] > 0: b_fore = mut_count_b[0]

        a_frac = mut_count_a / np.sum(df_wide[df_wide[subset_key] == 0][count_key].values)
        b_frac = mut_count_b / np.sum(df_wide[df_wide[subset_key] == 1][count_key].values)

        if mut_count_a.shape[0] == 0 or mut_count_b.shape[0] == 0: continue

        ratio = np.log2(b_frac / a_frac)[0]
        ratio = b_frac / a_frac

        odds[mut2idx[mut]] = ratio[0]

    if x is None: 
        x = odds
        continue
    if y is None: y = odds

f, ax = plt.subplots()

for mut in mut2idx:
    i = mut2idx[mut]
    ax.scatter(x[i], y[i], label=mut)

ax.legend()

f.savefig('o.png')
