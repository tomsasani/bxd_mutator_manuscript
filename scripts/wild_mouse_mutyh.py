import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns 
import argparse

p = argparse.ArgumentParser()
#p.add_argument("--strain_vars")
#p.add_argument("--out_a")
#p.add_argument("--out_b")

args = p.parse_args()



df = pd.read_csv("/Users/tomsasani/harrislab/bxd_mutator_ms/data/dumont.mutyh.csv", 
                    names=["interval", "strain", "gt", "ref", "alt"])

samps = pd.unique(df['strain'])
intervals = pd.unique(df['interval'])

mapping = pd.read_excel("/Users/tomsasani/harrislab/bxd_mutator_ms/data/wild_mouse_sample_list.xlsx")

mapping['strain'] = mapping['Sample'].astype(str).apply(lambda s: s.rstrip())

mapping = mapping.sort_values("Subspecies")

starts, ends = np.zeros(4), np.zeros(4)
idx = 0
prev_s = None
for i in np.arange(mapping.shape[0]):
    if prev_s is None: 
        prev_s = mapping['Subspecies'].values[i]
        continue
    if mapping['Subspecies'].values[i] != prev_s:
        ends[idx] = i
        idx += 1
        prev_s = mapping['Subspecies'].values[i]
        starts[idx] = i
        continue
    else: continue
ends[-1] = mapping.shape[0]
print (pd.unique(mapping['Subspecies']))
print (starts, ends)

s2sp = dict(zip(mapping["strain"], mapping['Subspecies']))

out_a = np.zeros((len(intervals), len(samps)))

smp2idx = dict(zip(mapping['strain'], range(len(mapping['strain']))))
int2idx = dict(zip(intervals, range(len(intervals))))

#uniq_sp = set([s.split('_')[0] for s in smp2idx])
uniq_sp = set([s2sp[s] for s in smp2idx])

for i,row in df.iterrows():
    i_i = int2idx[row['interval']]
    strain = row['strain']
    
    s_i = smp2idx[row['strain']]
    out_a[i_i, s_i] = row['gt']

# make heatmap in two parts to save space
sns.set_style('ticks')
f1, (ax1, ax2) = plt.subplots(1, 2, figsize=(14,2.5), sharex=False, gridspec_kw={'width_ratios': [1, 3]})
f2, (ax3, ax4) = plt.subplots(1, 2, figsize=(14,2.5), sharex=False, gridspec_kw={'width_ratios': [4, 1]})

cmap = ['white', 'gainsboro', 'lightsteelblue', 'royalblue']

interval2mut = dict(zip(intervals, ["p.Gln5Arg",
                                   "p.Arg24Cys",
                                   "p.Ser69Arg",
                                   "p.Thr312Pro",
                                   "p.Thr313Pro"]))

tick_colors_dict = dict(zip(uniq_sp, sns.color_palette('colorblind', len(uniq_sp) + 1)))

for ax, s, e, t in ((ax1, 0, 8, "$\it{M. spretus}$"), (ax2, 8, 38, '$\it{M. m. castaneus}$'), 
                    (ax3, 38, 140, '$\it{M. m. domesticus}$'), (ax4, 140, 154, '$\it{M. m. musculus}$')):

    cbar = False
    if ax == ax2 or ax == ax4: cbar = True
    sns.heatmap(out_a[:,s:e], ax=ax, cmap=cmap, linecolor='k', linewidth=1, cbar=cbar)

    if cbar:
        colorbar = ax.collections[0].colorbar
        colorbar.set_ticks([0.2, 0.8, 1.4, 2])
        colorbar.set_ticklabels(['unknown', 'homozgyous\nreference', 'heterozygous', 'homozgyous\nalternate'])

    #samps_reform = [s.split('_')[1] for s in samps][s:e]
    samps_reform = [s2sp[s] for s in samps][s:e]
    samps_reform = samps[s:e]
    samps_reform = mapping['strain'].values[s:e]
    ax.set_xticks(np.arange(len(samps_reform)) + 0.5)
    ax.set_xticklabels(samps_reform, rotation=90, fontsize=10)
    if ax == ax1 or ax == ax3:
        ax.set_yticks(np.arange(len(intervals)) + 0.5)
        ax.set_yticklabels([interval2mut[it] for it in intervals], rotation=0)
    else:
        ax.set_yticks([])

    #for i,ticklabel in enumerate(ax.get_xticklabels()):
    #    ticklabel.set_color(tick_colors_dict[samps[s:e][i].split('_')[0]])

    ax.set_title(t)
f1.tight_layout()
f2.tight_layout()
f1.savefig("/Users/tomsasani/harrislab/bxd_mutator_ms/o.png", bbox_inches='tight')
f2.savefig("/Users/tomsasani/harrislab/bxd_mutator_ms/oo.png", bbox_inches='tight')

#f2.savefig(args.out_b, bbox_inches='tight')

