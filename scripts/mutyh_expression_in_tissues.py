import pandas as pd 
import matplotlib.pyplot as plt 
import seaborn as sns
import glob
import argparse
import scipy.stats as ss

p = argparse.ArgumentParser()
p.add_argument("--tidy_spectra", required=True,
                    help="""file containing mutation spectra for \
                            each sample in tidy, long-form format""")
p.add_argument("--strain_metadata", required=True, 
                    help="""metadata about strains in Excel format.""")
p.add_argument("--out", required=True,
                    help="""name of output file""")
args = p.parse_args()

sns.set_style('ticks')

f, axarr = plt.subplots(2, 3, figsize=(8, 6))

# define X and Y axis coordinates for plot
y_coord = [0, 1, 2, 0, 1, 2]
x_coord = [0, 0, 0, 1, 1, 1]

cell_types = ["spleen", "hematopoietic_stem_cells", "liver",
              "retina", "kidney", "amygdala"]

for cell_type, x, y in zip(cell_types, x_coord, y_coord):

    fh = "data/gene_network_expression/{}_rnaseq.csv".format(cell_type)

    # read in expression data
    expr = pd.read_csv(fh, header=5, names=["gn_name", "expression", "a", "b"])

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

    sns.boxplot(x="haplotype_at_qtl", y="expression", ax=axarr[x, y], 
                    color='white', data=spectra, fliersize=0)
    sns.stripplot(x="haplotype_at_qtl", y="expression", ax=axarr[x, y], 
                        palette="colorblind", s=5, data=spectra)

    a_vals = spectra.query('haplotype_at_qtl == "D"')['expression'].values
    b_vals = spectra.query('haplotype_at_qtl == "B"')['expression'].values

    stat, p = ss.ttest_ind(a_vals, b_vals, equal_var=False)

    if y == 0:
        axarr[x, y].set_ylabel("Mutyh expression")
    else: axarr[x, y].set_ylabel(None)
    if x == 1:
        axarr[x, y].set_xlabel("Haplotype at rs52263933")
    else: axarr[x, y].set_xlabel(None)

    if x == 0:
        axarr[x, y].set_xticks([])
        axarr[x, y].set_xticklabels([])

    title = ' '.join(cell_type.split('_')).capitalize()
    if p > 0.01:
        title += '\n(Welch t-test p = {})'.format(round(p, 2))
    else:
        title += '\n(Welch t-test p = {:0.2e})'.format(p)
    axarr[x, y].set_title(title)

    sns.despine(ax=axarr[x, y])

f.tight_layout()

f.savefig(args.out, dpi=300, bbox_inches='tight')
