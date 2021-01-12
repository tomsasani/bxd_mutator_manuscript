import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns 

df = pd.read_csv('/Users/tomsasani/harrislab/bxd_mutator_ms/data/wild_mouse_mutyh.csv', names=["interval", "strain", "gt", "ref", "alt"])

samps = pd.unique(df['strain'])
intervals = pd.unique(df['interval'])

#sample_map = pd.read_csv("/Users/tomsasani/harrislab/bxd_mutator_ms/data/wild_mouse_sample_list.csv")
#sample_map['Sample'] = sample_map['Sample'].apply(lambda s: s.rstrip())

#samps = [s for s in samps if not s.startswith("SRR")]
#samps = sample_map['Sample'].values.astype(str)
#samps = 
out_a = np.zeros((len(intervals), len(samps)))

smp2idx = dict(zip(samps, range(len(samps))))
int2idx = dict(zip(intervals, range(len(intervals))))

uniq_sp = set([s.split('_')[0] for s in samps])
#uniq_sp = pd.unique(sample_map['Subspecies'])

sp_an = np.zeros((len(uniq_sp), len(intervals)))
sp_ac = np.zeros((len(uniq_sp), len(intervals)))
#smp2sp = dict(zip(sample_map['Sample'].values.astype(str), sample_map['Subspecies']))
sp2idx = dict(zip(uniq_sp, range(len(uniq_sp))))

for i,row in df.iterrows():
    i_i = int2idx[row['interval']]
    strain = row['strain']
    #if row['strain'] not in smp2idx: 
    #    strain = row['strain'] + u'\xa0'
     #   continue
    s_i = smp2idx[row['strain']]
    #species = smp2sp[row['strain']]
    sp_i = sp2idx[row['strain'].split('_')[0]]
    #sp_i = sp2idx[species]
    out_a[i_i, s_i] = row['gt']

    if row['gt'] == -1: continue
    sp_an[sp_i, i_i] += 2
    sp_ac[sp_i, i_i] += row['gt']

f, ax = plt.subplots( figsize=(14,2))

cmap = ['white', 'lightgrey', 'lightsteelblue', 'royalblue']

interval2mut = dict(zip(intervals, ["p.Gln5Arg",
                                   "p.Arg24Cys",
                                   "p.Ser69Arg",
                                   "p.Thr312Pro",
                                   "p.Thr313Pro"]))

sns.heatmap(out_a, ax=ax, cmap=cmap[1:], edgecolor='k', lw=0.1)

print (out_a)

samps_reform = ['_'.join(s.split('_')[:2]) for s in samps]
#samps_reform = samps
ax.set_xticks(np.arange(len(samps_reform)) + 0.5)
ax.set_yticks(np.arange(len(intervals)) + 0.5)
ax.set_yticklabels([interval2mut[it] for it in intervals], rotation=0)
ax.set_xticklabels(samps_reform, rotation=90)

f.savefig('heatmap.png', dpi=300, bbox_inches='tight')

f, ax = plt.subplots()

ind = np.arange(len(intervals))

colors = sns.color_palette('colorblind', len(uniq_sp))

sp_af = sp_ac / sp_an
print (sp_an)
print (sp_af)
width = 0.15
for i in np.arange(len(uniq_sp)):
    ax.bar(ind + (width * i), sp_af[i], width, color=colors[i], label=[s for s in sp2idx][i], edgecolor='k')

ax.legend(frameon=False)
ax.set_xticks(ind + 0.25)
ax.set_xticklabels([interval2mut[it] for it in intervals], rotation=30)


f.savefig('obar.png', dpi=300, bbox_inches='tight')