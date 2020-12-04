import matplotlib.pyplot as plt
import scipy.stats as ss
import numpy as np
import seaborn as sns
import pandas as pd
import numpy as np
import math
import argparse

def get_bootstrap_ci(cons, n_trials=100): 
    cons_trials = np.zeros((n_trials, cons.shape[0]))
    
    for t in range(n_trials):
        
        shuff_scores = np.random.choice(cons, cons.shape[0], replace=True)
        cons_trials[t, :] = shuff_scores   
    
    conf_ints = np.zeros((2, 99))

    # sort shuffled conservation scores by calculating difference
    # between empirical and true means
    cons_trials_diffs = np.mean(cons_trials, axis = 1) - np.mean(cons)
    
    sorted_idxs = cons_trials_diffs.argsort()
    sorted_vals = cons_trials[sorted_idxs]

    critval = int(0.99 * n_trials)
    lo_crit, hi_crit = n_trials - critval - 1, critval - 1

    lo_bound = sorted_vals[lo_crit]
    hi_bound = sorted_vals[hi_crit]

    conf_ints[0] = np.histogram(lo_bound, bins=99)[0]
    conf_ints[1] = np.histogram(hi_bound, bins=99)[0]
    
    conf_ints_fracs = conf_ints / np.sum(conf_ints, axis=1)[:,None]
    
    return conf_ints_fracs


p = argparse.ArgumentParser()
p.add_argument("--annotated_singletons")
p.add_argument("--annotated_common")
p.add_argument("--out")
args = p.parse_args()

singletons = pd.read_csv(args.annotated_singletons)
common = pd.read_csv(args.annotated_common)

singletons['v_type'] = 'singletons'
common['v_type'] = 'common'

combined = pd.concat([singletons, common])

f, ax = plt.subplots()

sns.set_style('ticks')

cols = sns.color_palette("colorblind", 2)

xticklabels = None

for i,lab in enumerate(pd.unique(combined['v_type'])):

    sub_df = combined[combined['v_type'] == lab]

    cons = sub_df[sub_df['phastCons'] != "None"]["phastCons"].values.astype(np.float64)
    
    y, x = np.histogram(cons, bins=99)
    y = y / np.sum(y)
    y = np.cumsum(y)

    ci = get_bootstrap_ci(cons) 

    ax.plot(x[1:], y, color=cols[i], label=lab.capitalize())
    ax.fill_between(x[1:], np.cumsum(ci[0]), np.cumsum(ci[1]), color=cols[i], alpha=0.25)

ax.legend(loc="lower right", title="Variant type", frameon=False)
ax.set_xlabel('Phastcons score (low {} high conservation)'.format(r'$\to$'))
ax.set_ylabel('Cumulative fraction of variants')
ax.set_xticks(np.arange(0, 1.1, 0.1))
#ax.set_xticklabels(xticklabels)
#ax.set_ylim(0, 1.05)

sns.despine(ax=ax, top=True, right=True)

f.savefig(args.out, bbox_inches='tight')
