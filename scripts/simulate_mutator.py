import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# read in tidy data
tidy_spectra = pd.read_csv("/Users/tomsasani/harrislab/bxd_mutator_ms/csv/tidy_mutation_spectra.csv")

tidy_spectra = tidy_spectra.query('estimate_type == "count"')

mut2idx = dict(zip(pd.unique(tidy_spectra['base_mut']), range(len(pd.unique(tidy_spectra['base_mut'])))))

muts_per_strain = np.zeros((len(pd.unique(tidy_spectra['bxd_strain_conv'])), len(mut2idx)), dtype=np.int64)

for s_i,s in enumerate(pd.unique(tidy_spectra['bxd_strain_conv'])):
    for mut in mut2idx:
        v = tidy_spectra[(tidy_spectra['bxd_strain_conv'] == s) & (tidy_spectra['base_mut'] == mut)]['estimate'].values
        if v.shape[0] == 0: continue
        muts_per_strain[s_i, mut2idx[mut]] += v[0]

n_samps = muts_per_strain.shape[0]

samp_sums = np.sum(muts_per_strain, axis=1)

mut_fracs = muts_per_strain / samp_sums[:,None]

d2_idxs = np.random.choice(np.arange(n_samps), size=int(n_samps / 2))

new_counts = np.zeros(mut_fracs.shape, dtype=np.float64)

for s_i in np.arange(mut_fracs.shape[0]):
    n_muts = samp_sums[s_i]
    new_sample = np.random.multinomial(n_muts, mut_fracs[s_i])
    new_counts[s_i] = new_sample

d2_counts = new_counts[d2_idxs]

mut2increase = "C>A"

idx2increase = mut2idx[mut2increase]

d2_counts[:,idx2increase] *= 1.5

d2 = pd.DataFrame(d2_counts, columns=list(mut2idx.keys()))

b6_counts = new_counts[~d2_idxs]

b6 = pd.DataFrame(b6_counts, columns=list(mut2idx.keys()))

d2['strain'] = 'D2'
b6['strain'] = 'B6'

combined = pd.concat([d2, b6])

combined['totals'] = combined.apply(lambda row: np.sum(row.values[:-1]), axis=1)

combined_tidy = combined.melt(id_vars=["strain", "totals"])
combined_tidy['frac'] = combined_tidy['value'] / combined_tidy['totals']

f, ax = plt.subplots()
sns.boxplot(x="variable", y="frac", hue="strain", data=combined_tidy, ax=ax)
f.savefig('o.png')



