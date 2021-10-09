import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import argparse
import scipy.stats as ss

plt.rcParams.update({'font.size': 16})

p = argparse.ArgumentParser()
p.add_argument("--tidy_spectra", required=True,
                    help="""file containing mutation spectra for \
                            each sample in tidy, long-form format""")
p.add_argument("--strain_metadata", required=True,
                    help="""metadata about strains in Excel format.""")
p.add_argument("--expression_data", required=True)
p.add_argument("--out", required=True,
                    help="""name of output file""")
args = p.parse_args()

sns.set_style('ticks')

f, ax = plt.subplots(figsize=(8, 6))

# read in expression data
expr = pd.read_csv(args.expression_data,
                   header=5,
                   names=[
                       "gn_name",
                       "expression",
                       "a",
                       "b",
                   ])
# read in tidy mutation spectra data for BXDs
spectra = pd.read_csv(args.tidy_spectra)
# read in BXD strain metadata
summary = pd.read_excel(args.strain_metadata)
# map the "original" BXD names to their gene network names
strain2genenet = dict(zip(summary['bam_name'], summary['GeneNetwork name']))
# add a column with gene network names to the tidy mutation spectra
spectra['gn_name'] = spectra['bxd_strain'].apply(lambda s: strain2genenet[s])
# only use one entry per BXD
spectra = spectra.query('estimate_type == "rate" & base_mut == "C>A"')

strain2expr = dict(zip(expr['gn_name'], expr['expression']))

# add a column to the spectra dataframe with the expression value of each strain
# assuming the strain has expression data from gene network
spectra['expression'] = spectra['gn_name'].apply(lambda s: strain2expr[s] if s in strain2expr else -1)
spectra = spectra.query('expression != -1')

sns.boxplot(
    x="haplotype_at_qtl",
    y="expression",
    ax=ax,
    color='white',
    data=spectra,
    fliersize=0,
)
sns.stripplot(
    x="haplotype_at_qtl",
    y="expression",
    ax=ax,
    palette="colorblind",
    s=5,
    data=spectra,
)

a_vals = spectra.query('haplotype_at_qtl == "D"')['expression'].values
b_vals = spectra.query('haplotype_at_qtl == "B"')['expression'].values

stat, p = ss.ttest_ind(a_vals, b_vals, equal_var=False)

ax.set_ylabel("Mutyh expression")
ax.set_xlabel("Haplotype at rs52263933")

title = ' '.join("Gastrointestinal RNA")
if p > 0.01:
    title += '\n(Welch t-test p = {})'.format(round(p, 2))
else:
    title += '\n(Welch t-test p = {:0.2e})'.format(p)
ax.set_title(title)

sns.despine(ax=ax)

f.tight_layout()
f.savefig(args.out, bbox_inches='tight')
