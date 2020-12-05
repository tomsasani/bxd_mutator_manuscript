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
p.add_argument("--annotated_singletons", required=True, 
                    help="""annotated singleton variants in extended BED format.""")
p.add_argument("--out_pref", required=True,
                    help="""prefix to prepend to output files. two files will \
                            be produced, one with information about individual
                            substitution types, and one with information about
                            overall mutation rates.""")
args = p.parse_args()

# read in the mutation file
singleton = pd.read_csv(args.annotated_singletons)

# columns to group the singletons by
group_cols = ["bxd_strain_conv", 
              "epoch", 
              "n_inbreeding_gens", 
              "haplotype_at_qtl",
              "n_intercross_gens", 
              "n_callable_bp", 
              "base_mut"]

# convert to grouped, wide-form dataframe
df_wide = singleton.groupby(group_cols).count().add_suffix('_count').reset_index()

# subset wide-form dataframe to relevant columns
group_cols.append('chrom_count')
df_wide = df_wide[group_cols]

# remove any samples without information about callable bp or generation numbers.
# samples that have been backcrossed at some point during their inbreeding
# are marked with a "NA" in that column and are removed here.
df_wide = df_wide[df_wide['n_inbreeding_gens'] != "NA"]
df_wide = df_wide[df_wide['n_callable_bp'] != "NA"]

df_wide.rename(columns={'chrom_count': 'count'}, inplace=True)

# map samples to their sums of total mutations, so we can calculate
# relative fractions of each mutation type
smp_sums = df_wide.groupby(['bxd_strain_conv', 
                            'n_callable_bp', 
                            'n_inbreeding_gens', 
                            'n_intercross_gens', 
                            'epoch']).sum().add_suffix('_sum').reset_index()

smp_sums.rename(columns = {'count_sum': 'total_muts'}, inplace=True)

# calculate overall mutation rates as the sum of singletons divided
# by the number of generations of inbreeding and the diploid number
# of base pairs that were "callable" in the strain
smp_sums['rate'] = smp_sums['total_muts'] / smp_sums['n_inbreeding_gens'] / smp_sums['n_callable_bp'] / 2
smp2sum = dict(zip(smp_sums['bxd_strain_conv'], smp_sums['total_muts']))

df_wide['total_muts'] = df_wide['bxd_strain_conv'].apply(lambda s: smp2sum[s])

# remove samples with very few singletons
df_wide = df_wide.query('total_muts >= 50')

# add a column to the dataframe with the relative fraction of each mutation type
df_wide['fraction'] = df_wide['count'] / df_wide['total_muts']

# and add a column with the rate of mutation for each mutation type
df_wide['rate'] = df_wide['count'] / df_wide['n_inbreeding_gens'] / df_wide['n_callable_bp'] / 2

# make a tidy dataframe such that there are three entries for every mutation type
# in each sample -- one with its count, one with its rate, and one with its fraction
df_tidy = df_wide.melt(id_vars=group_cols[:-1], var_name="estimate_type", value_name="estimate")

df_tidy = df_tidy.query('estimate_type != "total_muts"')

# output the tidy dataframe with info about individual mutation types
df_tidy.to_csv(args.out_pref + "tidy_mutation_spectra.csv", index=False)

# output the tidy dataframe with info about overall mutation rates
smp_sums[['bxd_strain_conv',
          'n_callable_bp',
          'n_inbreeding_gens',
          'n_intercross_gens',
          'epoch',
          'total_muts',
          'rate']].to_csv(args.out_pref + "tidy_mutation_rates.csv", index=False)

