import pandas as pd
import scipy.stats as ss
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
import numpy as np
import argparse
from figure_gen_utils import revcomp
from scipy.spatial.distance import cosine

plt.rc('font', size=16)

p = argparse.ArgumentParser()
p.add_argument("--cosmic_signature", required=True, 
                    help="""file containing mutation spectrum corresponding to SBS18.""")
p.add_argument("--ohno_data", required=True,
                    help="""supplementary file from ohno et al. (2014).""")
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

def revcomp(nuc):
    d = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

    return ''.join([d[n] for n in list(nuc)])[::-1] 

def convert_cosmic_mutation(row):
    context = row['Subtype']
    mutation = row['Type']

    changed_from = mutation.split('>')[0]
    changed_to = mutation.split('>')[1]

    e_kmer = context[0] + changed_to + context[-1]

    if changed_from == "T":
        context = revcomp(context)
        e_kmer = revcomp(e_kmer)

    return context + '>' + e_kmer

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
print (cosine_similarity)