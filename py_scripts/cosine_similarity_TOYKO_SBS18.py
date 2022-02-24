import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
from figure_gen_utils import convert_toyko_mutation, convert_cosmic_mutation
from scipy.spatial.distance import cosine

plt.rc('font', size=16)

p = argparse.ArgumentParser()
p.add_argument(
    "--cosmic_signature",
    required=True,
    help="""file containing mutation spectrum corresponding to SBS18.""",
)
p.add_argument(
    "--ohno_data",
    required=True,
    help="""supplementary file from ohno et al. (2014).""",
)
args = p.parse_args()

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

# read in the COSMIC signature
cosmic = pd.read_csv(args.cosmic_signature)

# convert COSMIC mutation notation to match mine
cosmic['kmer'] = cosmic.apply(convert_cosmic_mutation, axis=1)

merged_df = cosmic.merge(toyko_wide, on="kmer")

cosine_similarity = 1 - cosine(merged_df.SBS18_mm10.values, merged_df.frac)
