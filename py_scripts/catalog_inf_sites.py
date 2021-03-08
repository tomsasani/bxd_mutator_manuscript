from cyvcf2 import VCF
import gzip
from quicksect import IntervalTree
from collections import defaultdict
import csv
import argparse
import numpy as np
import os
import sys
from utils import get_good_idxs, make_interval_tree

def make_output(fh):
	fh = open(fh, 'w')
	return fh

def make_chr_output(fh1, fh2):
	"""
	if we don't output our results per chromosome, the
	resulting file of states x positions is HUGE. and if
	we separate out our per-chromosome results into two files
	(one with just the sites, and one with the GT x sample matrix),
	we can read in the first relatively small file, and read in the
	matrix with dtype=int8, saving loads of memory. then, it's
	relatively trivial to access the positions for each row of the matrix
	by row index.
	"""

	# break if the file already exists so we don't 
	# erroneously overwrite
	fh1, fh2 = make_output(fh1), make_output(fh2)
	
	# add headers
	print ('### states with respect to B6', file=fh2)
	print ('### 0 = B6, 2 = D2, 1 = UNK', file=fh2)

	return fh1, fh2

p = argparse.ArgumentParser()
p.add_argument("--vcf", required=True, 
					help="""path to BXD VCF file""")
p.add_argument("--chrom", required=True,
					help="""chromosome to limit analysis""")
p.add_argument("--fh_a", required=True,
					help="""name of output file one (positions of each informative site)""")
p.add_argument("--fh_b", required=True, 
					help="""name of output file two (BXD genotypes at each informative site)""")
p.add_argument("-exclude")
args = p.parse_args()

# --
# read in VCF
# --
vcf = VCF(args.vcf, gts012=True)

# --
# define some global variables
# --
B6 = '4512-JFI-0333_C57BL_6J_two_lanes_phased_possorted_bam'
D2 = '4512-JFI-0334_DBA_2J_three_lanes_phased_possorted_bam'

HOM_REF, HET, HOM_ALT, UNK = range(4)

# --
# make some dictionaries
# --
samples2use = vcf.samples
samples2use = [s for s in samples2use if 'BXD' in s]
samples2use.extend([B6, D2])

vcf.set_samples(samples2use)

smp2idx = dict(zip(vcf.samples, range(len(vcf.samples))))
idx2smp = {v:k for k,v in smp2idx.items()} 

# --
# iterate over VCF
# --
exclude = None
if args.exclude:
	exclude = make_interval_tree(args.exclude)

fh1, fh2 = make_chr_output(args.fh_a, args.fh_b)

vcf_h = vcf(args.chrom)

print (','.join([s for s in vcf.samples if 'BXD' in s]), file=fh2)

for i,v in enumerate(vcf_h):
	if exclude and len(exclude[v.CHROM].search(v.start, v.end)) > 0: continue
	
	# limit to passing SNVs
	if v.var_type != "snp": continue
	if v.FILTER not in (None, "PASS"): continue
	if len(v.ALT) > 1: continue
	if v.call_rate < 1: continue

	b6_idx, d2_idx = smp2idx[B6], smp2idx[D2]

	gts = v.gt_types

	# if this site isn't a clear difference between B6 and D2
	# we're not interested
	if UNK in (gts[b6_idx], gts[d2_idx]): continue
	if abs(gts[b6_idx] - gts[d2_idx]) != 2: continue

	b6_gt = gts[b6_idx]
	d2_gt = gts[d2_idx]

	gq = v.gt_quals

	# get a list of indices corresponding to sites
	# that are not UNK and pass GQ/DP filters
	non_unk = np.where(gts != UNK)
	hi_qual = np.where(gq >= 30) 

	good_sites = np.intersect1d(non_unk, hi_qual)

	if b6_idx not in good_sites: continue
	if d2_idx not in good_sites: continue

	out_array = np.ones(len(vcf.samples), dtype=np.int8)

	smp_idx = np.arange(len(vcf.samples))

	good_idx = np.intersect1d(smp_idx, good_sites)
	good_idx_bxd = good_idx[(good_idx != b6_idx) & (good_idx != d2_idx)]

	good_gt_bxd = gts[good_idx_bxd]

	match_b6 = good_gt_bxd == b6_gt
	match_d2 = good_gt_bxd == d2_gt

	match_b6_loc = good_idx_bxd[match_b6]
	match_d2_loc = good_idx_bxd[match_d2]

	out_array[match_b6_loc] = 0
	out_array[match_d2_loc] = 2
	
	out_array = list(map(str, list(out_array)))

	print (str(v.POS), file=fh1)
	print (','.join(out_array), file=fh2)


