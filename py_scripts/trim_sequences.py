from skbio.alignment import TabularMSA
from skbio import DNA
from io import StringIO
import argparse
import numpy as np
from collections import Counter

p = argparse.ArgumentParser()
p.add_argument("--msa")
p.add_argument("-gap_frac", default=0.5, type=float)
p.add_argument("-only_plot_mutyh", action="store_true")
args = p.parse_args()

align = TabularMSA.read(args.msa, constructor=DNA, format="fasta" )

# we'll only plot the 6 amino acids of interest
aa_positions_to_keep = [5, 24, 69, 153, 312, 313]
nuc_positions_to_keep = []
for aa in aa_positions_to_keep:
    for nuc_i in range(3):
        nuc_position = (aa * 3) - 3 + nuc_i
        nuc_positions_to_keep.append(nuc_position)

bl6_seq_idx = 1
bl6_seq = str(align[bl6_seq_idx])

if not args.only_plot_mutyh:
    nuc_positions_to_keep = range(len(bl6_seq))

good_nucs = []

good_nuc_counter = 0
for nuc_i, nuc in enumerate(bl6_seq):
    if nuc == "-": continue
    else: 
        if good_nuc_counter in nuc_positions_to_keep:
            good_nucs.append(nuc_i)
        good_nuc_counter += 1

n_nucs = align.shape.position
n_seqs = align.shape.sequence

seq_to_keep = np.zeros(n_nucs, dtype=np.int8)

for nuc_i in np.arange(n_nucs):
    nucs_at_site = str(align[:,nuc_i])

    n_gaps = Counter(nucs_at_site)['-']

    if nuc_i not in good_nucs: continue
    if (n_seqs - n_gaps) / n_seqs >= args.gap_frac: seq_to_keep[nuc_i] = 1

with StringIO() as fh:
    print(align[:,seq_to_keep.astype(bool)].write(fh).getvalue())
