import matplotlib.pyplot as plt
import pandas as pd
import csv
import scipy.stats as ss
from collections import defaultdict, Counter
from quicksect import IntervalTree
import numpy as np
import seaborn as sns
import glob
import sys
import argparse

p = argparse.ArgumentParser()
p.add_argument("--annotated_singletons")
p.add_argument("--outdir")
args = p.parse_args()

# ---
# read in the mutation file
# ---

singleton = pd.read_csv(args.annotated_singletons)

untidy_cols = ["bxd_strain_conv", "epoch", "n_inbreeding_gens", "haplotype_at_qtl",
             "n_intercross_gens", "n_callable_bp", "base_mut"]

# convert to untidy dataframe
singleton_untidy = singleton.groupby(untidy_cols).count().add_suffix('_count').reset_index()

# subset untidy dataframe to relevant columns
untidy_cols.append('chrom_count')

singleton_untidy = singleton_untidy[untidy_cols]

# remove any samples without information about callable bp or generation numbers
singleton_untidy = singleton_untidy[singleton_untidy['n_inbreeding_gens'] != "NA"]
singleton_untidy = singleton_untidy[singleton_untidy['n_callable_bp'] != -1]

singleton_untidy.rename(columns={'chrom_count': 'count'}, inplace=True)

# map samples to their sums of total mutations, so we can calculate
# fractions of each mutation type
smp_sums = singleton_untidy.groupby(['bxd_strain_conv', 'n_callable_bp', 'n_inbreeding_gens', 
                                     'n_intercross_gens', 'epoch']).sum().add_suffix('_sum').reset_index()

smp_sums.rename(columns = {'count_sum': 'total_mutations'}, inplace=True)
smp_sums['rate'] = smp_sums['total_mutations'] / smp_sums['n_inbreeding_gens'] / smp_sums['n_callable_bp'] / 2
smp2sum = dict(zip(smp_sums['bxd_strain_conv'], smp_sums['total_mutations']))

smp_sums[['bxd_strain_conv', 'total_mutations', 'n_callable_bp', 'n_inbreeding_gens', 
          'n_intercross_gens', 'epoch', 'rate']].to_csv(args.outdir + "tidy_mutation_rates.csv", index=False)

singleton_untidy['total_mutations'] = singleton_untidy['bxd_strain_conv'].apply(lambda s: smp2sum[s])
singleton_untidy['fraction'] = singleton_untidy['count'] / singleton_untidy['total_mutations']


singleton_untidy['rate'] = singleton_untidy['count'] / singleton_untidy['n_inbreeding_gens'] / singleton_untidy['n_callable_bp'] / 2

singleton_untidy.to_csv(args.outdir + "tidy_mutation_fractions.csv", index=False)
