import argparse
import itertools
import numpy as np
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
import matplotlib

font = {'size'   : 16}

matplotlib.rc('font', **font)

def find_groups(a):
    """
    function to get the indexes of shared "groups"
    in an arbitrary list
    """
    groups = []
    cur_val, last_idx = 0, 0
    for idx,val in enumerate(a):
        if idx == 0:
            cur_val = val
            continue
        if val != cur_val or idx == len(a) - 1:
            if idx == len(a) - 1:
                groups.append((last_idx, idx + 1, cur_val))
            else:
                groups.append((last_idx, idx, cur_val))
            last_idx = idx
            cur_val = val
        else:
            cur_val = val
    return groups

p = argparse.ArgumentParser()
p.add_argument("--annotated_vars", required=True,
                help="""annotated variants in extended BED format""")
p.add_argument("--out", required=True,
                help="""name of output file""")
args = p.parse_args()

variants = pd.read_csv(args.annotated_vars)

# get a mapping of each strain to a corresponding index
smps = list(pd.unique(variants['bxd_strain_conv']))
smp2idx = dict(zip(smps, range(len(smps))))
idx2smp = {v:k for k,v in smp2idx.items()}

# get a mapping of each strain to its epoch of origin
smp2epoch = defaultdict()
for smp in smps:
    if smp in smp2epoch: continue
    epoch = pd.unique(variants[variants['bxd_strain_conv'] == smp]['epoch'])[0]
    smp2epoch[smp] = epoch

epochs = [smp2epoch[s] for s in smps]

# make a single column representing the VCF site
variants['site'] = variants['chrom'] + ':' + variants['start'].astype(str) + '-' + variants['end'].astype(str)

n_samps = len(set(variants['bxd_strain_conv']))
n_sites = len(set(variants['site']))

# generate a sample by sample matrix
smp_by_smp = np.zeros((n_samps, n_samps), dtype=np.int32)

# also generate a site by sample matrix
site_by_sample = np.zeros((n_sites, n_samps), dtype=np.int32)

# loop over every variant
for site_idx, site in enumerate(set(variants['site'])):

    # subset the dataframe to only include variants observed
    # in multiple samples at the same site
    shared_vars = variants[variants['site'] == site]

    # get a list of samples that share this variant
    smps_sharing = shared_vars['bxd_strain_conv'].values

    # get a list of haplotypes in the samples that share this variant
    haps_shared = shared_vars['haplotype'].values
    # and require all samples have the same haplotype
    if not len(set(haps_shared)) == 1: continue

    # increment the site x sample matrix by 1 for each
    # sample at this site
    for smp in smps_sharing:
        s_idx = smp2idx[smp]
        site_by_sample[site_idx, s_idx] += 1

    # for every pair of samples that share this variant,
    # increment the sample x sample matrix
    for s1, s2 in itertools.combinations(smps_sharing, 2):
        s1_idx, s2_idx = smp2idx[s1], smp2idx[s2]
        # increment both directions of the pair for symmetry
        smp_by_smp[s1_idx, s2_idx] += 1
        smp_by_smp[s2_idx, s1_idx] += 1

# construct a new pairwise matrix that describes the number of variants
# at which two samples *could have* shared a variant. in other words,
# the union of two sets: the sites where sample A shares a variant with
# any other sample, and the sites where sample B shares a variant with 
# any other sample. 
sharing_index = np.zeros((n_samps, n_samps))

for s1_idx in np.arange(n_samps):
    for s2_idx in np.arange(n_samps):
        if s1_idx == s2_idx: sharing_index[s1_idx, s2_idx] = 1

        # get a list of sites where sample 1 shares a variant
        # with any other sample. do the same for sample 2
        s1_sites = site_by_sample[:,s1_idx]
        s2_sites = site_by_sample[:,s2_idx]

        s1_nonzero = np.where(s1_sites > 0)[0]
        s2_nonzero = np.where(s2_sites > 0)[0]

        shared_nonzero = np.union1d(s1_nonzero, s2_nonzero)

        sharing_index[s1_idx, s2_idx] = shared_nonzero.shape[0]

# convert the pairwise numbers of shared sites to fractions by dividing each 
# pair's number of shared variants by the total number of unique variants that
# those two samples shared with all other samples
#pairwise_sharing_frac = np.nan_to_num(smp_by_smp / sharing_index)
pairwise_sharing_frac = smp_by_smp / sharing_index

# sort strains by the epoch from which they were derived
sorted_epoch_idxs = np.argsort(epochs)

# then, sort the pairwise sharing matrix by epoch so that we can
# plot samples of similar epochs next to each other
pairwise_sharing_sorted = pairwise_sharing_frac[sorted_epoch_idxs][:,sorted_epoch_idxs]

f, ax = plt.subplots(figsize=(12, 9))

# plot the heatmap of pairwise sharing
sns.heatmap(pairwise_sharing_sorted, ax=ax, vmin=0, vmax=0.5)#, cmap="bwr")

# find the sample indexes that correspond to the ends
# of epoch groups
epoch_idxs = find_groups(np.array(epochs)[sorted_epoch_idxs])

ax_ticks, ax_labels = [], []

# use coordinates that demarcate starts and ends
# of epoch groups to draw bounding boxes around heatmap
# cells that correspond to each epoch
for xy in epoch_idxs:
    x, y, e = xy
    size = y - x
    
    rect = patches.Rectangle((x,x), size, size, linewidth=1, 
                                edgecolor='w', facecolor='none', 
                                label='Epoch {}'.format(e))

    ax.add_patch(rect)
    ax_ticks.append((y + x) / 2)
    ax_labels.append(e)

# figure formatting
ax.set_yticks(ax_ticks)
ax.set_yticklabels(ax_labels)

ax.set_xticks(ax_ticks)
ax.set_xticklabels(ax_labels, rotation=0)

ax.set_xlabel("Epoch")
ax.set_ylabel("Epoch")

ax.set_title("Degree of BXD-private mutation sharing between samples")

f.savefig(args.out, dpi=300, bbox_inches='tight')
