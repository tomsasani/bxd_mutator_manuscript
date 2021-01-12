import argparse
import scipy.stats as ss
import itertools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def revcomp(kmer):

    d = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

    a, b = kmer.split('>')

    a_new = ''.join([d[n] for n in list(a)][::-1])
    b_new = ''.join([d[n] for n in list(b)][::-1])

    return '{}>{}'.format(a_new, b_new)

p = argparse.ArgumentParser()
p.add_argument("--dumont_xls")
args = p.parse_args()

dumont = pd.read_excel(args.dumont_xls, 
                       sheet_name="TableS4",
                       header=2)

dumont = dumont[dumont['FocalStrain'].isin(["DBA_2J", "C57BL_6NJ"])]

dumont['kmer_old'] = dumont.apply(lambda row: "{}{}{}>{}{}{}".format(row['5prime_flank'],
                                                                 row['ancestral'],
                                                                 row['3prime_flank'],
                                                                 row['5prime_flank'],
                                                                 row['derived'],
                                                                 row['3prime_flank']), axis=1)
#dumont_wide = dumont.groupby(['FocalStrain', 'kmer']).count().add_suffix('_count').reset_index()

dumont = dumont[['FocalStrain', 'kmer_old']]


dumont['kmer'] = dumont['kmer_old'].apply(lambda k: revcomp(k) if k[1] in ('T', 'G') else k)

dumont['value'] = 1

dumont.rename(columns={'FocalStrain':'strain'}, inplace=True)

dumont.to_csv("csv/dumont_singletons.csv", index=False)
