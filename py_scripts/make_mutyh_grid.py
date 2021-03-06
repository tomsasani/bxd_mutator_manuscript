import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

p = argparse.ArgumentParser()
p.add_argument(
    "--strain_vars",
    required=True,
    help="""genotypes of strains at each of the five
                        Mutyh missense mutations""",
)
p.add_argument(
    "--out",
    required=True,
    help="""name of output plot""",
)
p.add_argument(
    "-is_wild",
    action="store_true",
    help="""specify if samples are wild mice""",
)

args = p.parse_args()

wild = args.is_wild

if wild: plt.rc("font", size=32)
else: plt.rc("font", size=24)

df = pd.read_csv(args.strain_vars,
                 names=["interval", "strain", "gt", "ref", "alt"])

# remove strains that don't have associated singleton
# data from Dumont (2019)
not_in_dumont = [
    "CAST_EiJ",
    "LEWES_EiJ",
    "MOLF_EiJ",
    "PWK_PhJ",
    "SPRET_EiJ",
    "WSB_EiJ",
    "ZALENDE_EiJ",
]

df = df[~df['strain'].isin(not_in_dumont)]

samps = pd.unique(df['strain'])
intervals = pd.unique(df['interval'])

if wild:
    species = [s.split('_')[0] for s in samps]
else:
    species = ["Mouse Genome Project strains" for s in samps]

# keep track of unique species in dataset
uniq_sp = []
for s in species:
    if s in uniq_sp: continue
    uniq_sp.append(s)

# store the indices that correspond to the
# starts and ends of each species group in the
# dataframe of samples * genotypes
starts = np.zeros(len(uniq_sp), dtype=np.int8)
ends = np.zeros(len(uniq_sp), dtype=np.int8)
ends[-1] = len(species)

idx = 0
prev_s = None
for i, s in enumerate(species):
    # increment start and end index
    # when we see a new species in the list
    if prev_s is None:
        prev_s = s
        continue
    if s != prev_s:
        ends[idx] = i
        idx += 1
        prev_s = s
        starts[idx] = i
        continue
    else:
        continue

# create array of size n_intervals * n_samps
# to store genotypes for each strain
out_a = np.zeros((len(intervals), len(samps)))

# map samples and intervals to corresponding indices
# for easy array assignment
smp2idx = dict(zip(samps, range(len(samps))))
int2idx = dict(zip(intervals, range(len(intervals))))

# map intervals to actual a.a. changes
interval2mut = dict(
    zip(intervals, [
        "p.Gln5Arg",
        "p.Arg24Cys",
        "p.Ser69Arg",
        "p.Thr312Pro",
        "p.Ser313Pro",
    ]), )

# convert dataframe into `out_a` array of
# samples by genotypes
for i, row in df.iterrows():
    i_i = int2idx[row['interval']]
    strain = row['strain']

    s_i = smp2idx[row['strain']]
    out_a[i_i, s_i] = row['gt']

sns.set_style('ticks')

sp2fullname = {
    "Mmc": "castaneus",
    "Mmd": "domesticus",
    "Ms": "spretus",
    "Mmm": "musculus",
}

# define colormap for heatmap grid
cmap = ['gainsboro', 'lightsteelblue', 'royalblue']
p_i = 0

if wild:
    f, axarr = plt.subplots(
        1,
        4,
        sharey=True,
        figsize=(20, 4.5),
        gridspec_kw={'width_ratios': [4, 6.5, 0.75, 0.75]},
    )
else:
    f, axarr = plt.subplots(1, 1, figsize=(20, 3))

# define colors for MGP groups
colors = sns.color_palette('colorblind', 3)

# for every "group" of samples in a subspecies, generate a heatmap
for s, e, sp in zip(starts, ends, uniq_sp):
    # get genotypes of samples in the subspecies
    values = out_a[:, s:e]
    # get the unique haplotypes observed in those samples,
    # as well as the haplotype frequencies
    unique_haps, counts = np.unique(values, axis=1, return_counts=True)
    # sort haplotypes in order of frequency
    sorted_idxs = np.argsort(counts)[::-1]

    # set plotting variables, which will be different
    # if we're making grids for the wild vs. MGP mice
    # set plotting variables, which will be different
    # if we're making grids for the wild vs. MGP mice
    grid_vals = values
    hap_sums = np.sum(grid_vals, axis=0)
    sorted_haps = np.argsort(hap_sums)[::-1]
    grid_vals = grid_vals[:, sorted_haps]

    ax2use = axarr
    cmap = ['gainsboro', 'royalblue']

    lw = 2
    xticks = np.arange(values.shape[1]) + 0.5
    xlabs = [l.replace('_', '/') for l in samps[sorted_haps]]
    xlab_rotation = 90
    if wild:
        cmap = ['gainsboro', 'lightsteelblue', 'royalblue']

        ax2use = axarr[p_i]
        grid_vals = unique_haps[:, sorted_idxs]
        lw = 3
        xticks = np.arange(unique_haps.shape[1]) + 0.5
        xlabs = counts[sorted_idxs]
        xlab_rotation = 0
        title_rotation = 20

    sns.heatmap(grid_vals,
                ax=ax2use,
                cmap=cmap,
                linecolor='k',
                linewidth=lw,
                cbar=False)

    ax2use.set_xticks(xticks)
    ax2use.set_xticklabels(xlabs, rotation=xlab_rotation)

    if not wild:
        for tick_i, ticklabel in enumerate(ax2use.get_xticklabels()):
            if np.sum(grid_vals[:, tick_i]) == 10:
                ticklabel.set_color(colors[0])
            elif np.sum(grid_vals[:, tick_i]) == 6:
                ticklabel.set_color(colors[1])
            elif np.sum(grid_vals[:, tick_i]) == 0:
                ticklabel.set_color(colors[2])
            else:
                continue

    p_i += 1

# separate formatting for wild vs. MGP plots
if wild:
    axarr[1].set_xlabel("Number of samples with haplotype")
    axarr[0].set_yticks(np.arange(len(intervals)) + 0.5)
    axarr[0].set_yticklabels([interval2mut[it] for it in intervals],
                             rotation=0)
else:
    axarr.set_yticks(np.arange(len(intervals)) + 0.5)
    axarr.set_yticklabels([interval2mut[it] for it in intervals], rotation=0)

f.savefig(args.out, bbox_inches='tight')
