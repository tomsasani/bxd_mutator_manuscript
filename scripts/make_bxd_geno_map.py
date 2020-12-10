import pandas as pd
import numpy as np
import argparse

def convert_bxd_name(name: str) -> str:
    """
    depending on the VCF file we're using we may
    have to convert BXD strain names for portability
    """
    bxd_line_name = '_'.join(name.split('_phased')[0].split('_')[1:])
    bxd_line_num = name.split('_phased')[0].split('_')[0].split('-')[-1]

    bxd_line_new = bxd_line_name + '_' + bxd_line_num

    return bxd_line_new

p = argparse.ArgumentParser()
p.add_argument("--geno_file")
p.add_argument("--output")
args = p.parse_args()

# read in the genotypes file
geno = pd.read_csv(args.geno_file)

geno = geno.rename(columns={'Unnamed: 0': 'marker'})

# get a new header with converted BXD names in order to
# match names in singleton mutation dataframes
header = list(geno)
header_conv = [convert_bxd_name(n) for n in header[1:]]

# access the genotypes at each site
gts = geno.values[:,1:]

def gt_to_allele(gt):
    d = {0: 'B', 1: 'H', 2: 'D', -1:'NA'}
    return d[gt]

# convert genotypes to haplotypes
vf = np.vectorize(gt_to_allele)
new_gts = vf(gts)

# make a new dataframe with haplotypes
new_geno = pd.DataFrame(new_gts, columns=header_conv)
new_geno['marker'] = geno['marker']

# get indices of sites where all haplotypes are UNK
n_samples_unk_per_site = np.sum((new_geno == "NA").values, axis=1)
sites_all_unk = np.where(n_samples_unk_per_site == new_geno.shape[1] - 1)[0]
sites_all_not_unk = np.where(n_samples_unk_per_site < new_geno.shape[1] - 1)[0]

# subset to sites where at least one sample has a called haplotype
hq_markers = new_geno.iloc[sites_all_not_unk]

# make a new header
reordered_header = ['marker']
reordered_header.extend([h for h in header_conv if h != "NA"])

hq_markers_reordered = hq_markers[reordered_header]

hq_markers_reordered.to_csv(args.output, index=False)
