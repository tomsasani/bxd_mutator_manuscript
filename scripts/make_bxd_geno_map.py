import pandas as pd
import numpy as np
import argparse
from utils import convert_bxd_name

def gt_to_allele(gt: int) -> int:
    """
    convert integer genotypes to haplotypes
    """
    d = {0: 'B', 1: 'H', 2: 'D', -1:'NA'}
    return d[gt]

p = argparse.ArgumentParser()
p.add_argument("--geno_file", required=True, 
                    help="""file containing genotypes at roughly 7,000 markers for each \
                            BXD strain""")
p.add_argument("--output", required=True,
                    help="""name of output file""")
args = p.parse_args()

# read in the genotypes file
geno = pd.read_csv(args.geno_file)

geno = geno.rename(columns={'Unnamed: 0': 'marker'})

# remove founder samples from the output file
founder_samps = ['4512-JFI-0333_C57BL_6J_phased_possorted_bam',
                 '4512-JFI-0333_C57BL_6J_two_lanes_phased_possorted_bam',
                 '4512-JFI-0334_DBA_2J_phased_possorted_bam',
                 '4512-JFI-0334_DBA_2J_three_lanes_phased_possorted_bam']

geno = geno.drop(columns=founder_samps)

# get a new header with converted BXD names in order to
# match names in singleton mutation dataframes
header = list(geno)
header_conv = [convert_bxd_name(n) for n in header[1:]]

# access the genotypes at each site
gts = geno.values[:,1:]

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