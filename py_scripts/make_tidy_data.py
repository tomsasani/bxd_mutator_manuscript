import pandas as pd
import scipy.stats as ss
from collections import defaultdict
import numpy as np
import argparse
from figure_gen_utils import calc_new_gens, clr

p = argparse.ArgumentParser()
p.add_argument(
    "--annotated_singletons",
    required=True,
    help="""annotated singleton variants in extended BED format.""",
)
args = p.parse_args()

# read in the mutation file
singleton = pd.read_csv(args.annotated_singletons)

smp2hap = dict(zip(singleton.bxd_strain_conv, singleton.haplotype_at_qtl))

# columns to group the singletons by
group_cols = [
    "bxd_strain_conv",
    "bxd_strain",
    "epoch",
    "n_inbreeding_gens",
    "haplotype_at_qtl",
    "n_intercross_gens",
    "n_callable_bp",
    "base_mut",
]

# convert to grouped, wide-form dataframe
df_wide = singleton.groupby(group_cols).count().add_suffix(
    '_count').reset_index()

# subset wide-form dataframe to relevant columns
group_cols.append('chrom_count')
df_wide = df_wide[group_cols]

df_wide.rename(columns={'chrom_count': 'count'}, inplace=True)

# map samples to their sums of total mutations, so we can calculate
# relative fractions of each mutation type
smp_sums = df_wide.groupby([
    'bxd_strain_conv',
    'n_callable_bp',
    'bxd_strain',
    'n_inbreeding_gens',
    'n_intercross_gens',
    'epoch',
]).sum().add_suffix('_sum').reset_index()

smp_sums.rename(columns={'count_sum': 'total_muts'}, inplace=True)

print("Total of {} mutations in {} strains".format(
    np.sum(smp_sums['total_muts']),
    len(pd.unique(smp_sums['bxd_strain_conv']))))

# calculate overall mutation rates as the sum of singletons divided
# by the number of generations of inbreeding and the diploid number
# of base pairs that were "callable" in the strain
smp_sums['l_n'] = smp_sums['n_inbreeding_gens'].apply(
    lambda n: calc_new_gens(n))
smp_sums['rate'] = smp_sums['total_muts'] / (smp_sums['l_n'] *
                                             smp_sums['n_callable_bp'])
smp2sum = dict(zip(smp_sums['bxd_strain_conv'], smp_sums['total_muts']))

mean_mutation_rate = np.mean(smp_sums['rate'])
stderr_mutation_rate = ss.sem(smp_sums['rate'])

print("Mean mutation rate across strains is {} (95% CI = {} - {})".format(
    mean_mutation_rate,
    mean_mutation_rate - (2 * stderr_mutation_rate),
    mean_mutation_rate + (2 * stderr_mutation_rate),
))
d_smps = [s for s in smp2hap if smp2hap[s] == "D"]
b_smps = [s for s in smp2hap if smp2hap[s] == "B"]

mean_d_mutation_rate = np.mean(
    smp_sums[smp_sums['bxd_strain_conv'].isin(d_smps)].rate)
mean_b_mutation_rate = np.mean(
    smp_sums[smp_sums['bxd_strain_conv'].isin(b_smps)].rate)

print("Strains with D haplotypes accumulate mutations at {} the rate".format(
    mean_d_mutation_rate / mean_b_mutation_rate))

df_wide['total_muts'] = df_wide['bxd_strain_conv'].apply(lambda s: smp2sum[s])

# add a column to the dataframe with the relative fraction of each mutation type
df_wide['fraction'] = df_wide['count'] / df_wide['total_muts']

samps = pd.unique(df_wide['bxd_strain_conv'])
muts = pd.unique(df_wide['base_mut'])

smp_by_frac = np.zeros((len(samps), len(muts)))

for s_i, s in enumerate(samps):
    for mut_i, m in enumerate(muts):
        val = df_wide[(df_wide['bxd_strain_conv'] == s) & \
                (df_wide['base_mut'] == m)]['fraction'].values
        if val.shape[0] == 0: continue
        smp_by_frac[s_i, mut_i] = val[0]


sbf_clr = clr(smp_by_frac)

smp2mut2clr = defaultdict(lambda: defaultdict(float))
for s_i, s in enumerate(samps):
    for mut_i, m in enumerate(muts):
        smp2mut2clr[s][m] = sbf_clr[s_i, mut_i]

# and add a column with the rate of mutation for each mutation type
df_wide['l_n'] = df_wide['n_inbreeding_gens'].apply(lambda n: calc_new_gens(n))
df_wide['rate'] = df_wide['count'] / (df_wide['l_n'] *
                                      df_wide['n_callable_bp'])

df_wide['clr_fraction'] = df_wide.apply(
    lambda row: smp2mut2clr[row['bxd_strain_conv']][row['base_mut']], axis=1)

df_wide.fillna(value=0, inplace=True)

df_wide.replace(to_replace=np.inf, value=0, inplace=True)

# make a tidy dataframe such that there are three entries for every mutation type
# in each sample -- one with its count, one with its rate, and one with its fraction
df_tidy = df_wide.melt(
    id_vars=group_cols[:-1],
    var_name="estimate_type",
    value_name="estimate",
)

df_tidy = df_tidy.query('estimate_type != "total_muts"')

print(
    df_tidy.query('estimate_type == "count"').groupby(
        'base_mut').sum().reset_index())

for hap in ["D", "B"]:
    sub_df = df_tidy.query(
        "haplotype_at_qtl == @hap & base_mut == 'C>A' & estimate_type == 'rate'"
    )

    mean_mutation_rate = np.mean(sub_df['estimate'])
    stderr_mutation_rate = ss.sem(sub_df['estimate'])

    print(
        "Mean C>A mutation rate across strains with {} haps is {} (95% CI = {} - {})"
        .format(
            hap,
            mean_mutation_rate,
            mean_mutation_rate - (2 * stderr_mutation_rate),
            mean_mutation_rate + (2 * stderr_mutation_rate),
        ))

# output the tidy dataframe with info about individual mutation types
df_tidy.to_csv("csv/tidy_mutation_spectra.csv", index=False)

# output the tidy dataframe with info about overall mutation rates
smp_sums[[
    'bxd_strain_conv',
    'n_callable_bp',
    'n_inbreeding_gens',
    'bxd_strain',
    'l_n',
    'n_intercross_gens',
    'epoch',
    'total_muts',
    'rate',
]].to_csv("csv/tidy_mutation_rates.csv", index=False)
