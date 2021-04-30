import tabix
import os
import sys
import numpy as np
import gzip
import csv
import argparse
from cyvcf2 import VCF
from mutyper.ancestor import Ancestor
import doctest
from singleton_calling_utils import *

def run(args):

	# -----------
	# read in VCF
	# -----------
	vcf = VCF(args.vcf, gts012=True)

	# ------------
	# define some global variables that will come in handy
	# ------------
	HOM_REF, HET, HOM_ALT, UNK = range(4)

	B6 = '4512-JFI-0333_C57BL_6J_two_lanes_phased_possorted_bam'
	D2 = '4512-JFI-0334_DBA_2J_three_lanes_phased_possorted_bam'

	# -----------
	# generate a few necessary dictionaries
	# -----------

	# map the reference chromosome to a  
	# mutyper Ancestor() object
	ancestor = Ancestor(args.ref, k=args.nmer, sequence_always_upper=True)

	# subset samples to remove one of each pair of isogenic BXDs
	isogenic_smps = ["4512-JFI-0348_BXD24_TyJ_Cep290_J_phased_possorted_bam", # iso with BXD24
					  "4512-JFI-0344_BXD29_Ty_phased_possorted_bam", # iso with BXD29
					  "4512-JFI-0355_BXD152_phased_possorted_bam", # iso with BXD 155
					  "4512-JFI-0482_BXD087_RwwJ_phased_possorted_bam", # iso with BXD 194
					  "4512-JFI-0382_BXD048a_RwwJ_phased_possorted_bam", # iso with BXD48
					  "4512-JFI-0388_BXD65b_RwwJ_phased_possorted_bam", # iso with BXD65
					  "4512-JFI-0387_BXD65a_RwwJ_phased_possorted_bam", # iso with BXD65
					  "4512-JFI-0440_BXD073b_RwwJ_phased_possorted_bam", # iso with BXD73
					  "4512-JFI-0439_BXD73a_RwwJ_phased_possorted_bam" # iso with BXD73
					  ]
	samples2use = []
	for smp in vcf.samples:
		if "BXD" not in smp: continue
		if smp in isogenic_smps: continue
		samples2use.append(smp)

	samples2use.extend([B6, D2])
	vcf.set_samples(samples2use)

	samples_in_vcf = list(vcf.samples)

	# map sample indices to sample IDs and vice versa
	idx2smp = dict(zip(range(len(samples_in_vcf)), samples_in_vcf))
	smp2idx = dict(zip(samples_in_vcf, range(len(samples_in_vcf))))

	# ------------
	# generate various output files
	# ------------
	outcsvfh = open(args.out, 'w')

	csvheader = ['chrom', 'start', 'end', 'bxd_strain', 'kmer', 
				 'haplotype', 'gt', 'dp', 'ab', 'phastCons']

	print (','.join(csvheader), file=outcsvfh)

	# ------------
	# iterate over the VCF
	# ------------
	vcf_h = vcf(args.chrom)

	# generate an exclude file if it's provided
	exclude = None
	if args.exclude:
		exclude = make_interval_tree(args.exclude)
	
	# read in phastCons file
	tb = tabix.open(args.pcons)

	for i,v in enumerate(vcf_h):
		# keep track of reasons for each variant getting filtered out
		# check if the current variant is in the exclude file
		if exclude and len(exclude[v.CHROM].search(v.start, v.end)) > 0: continue

		if v.FILTER not in (None, "PASS"): continue
		if '*' in v.ALT: continue
		# limit to SNVs
		if v.var_type != "snp": continue
		if len(v.REF) > 1: continue

		for alt_idx, alt in enumerate(v.ALT):
			if len(alt) > 1: continue

			ref_allele, alt_allele = v.REF, alt

			# if the REF allele is longer than one nucleotide, but
			# the REF and ALT alleles are the same length, we almost
			# certainly need to normalize to become an SNV
			if len(v.REF) > 1 and (len(v.REF) == len(alt)):
				ref_allele, alt_allele = normalize_var(v.REF, alt)
				if 'N' in (ref_allele, alt_allele): continue
			
			# define the alternate genotype at this site.
			# if this is a multiallelic site, the genotype (integer)
			# corresponding to the first allele is 1, the genotype
			# corresponding to the second allele is 2, and so on.
			alt_gt = alt_idx + 1

			# access sample genotypes, excluding boolean flag indicating
			# whether sample is phased
			gts = np.array(v.genotypes)[:,:-1]

			# reformat genotypes
			gts_reformatted = reformat_genotypes(gts, alt_gt=alt_gt)
			

			# the index of the reference allele is always 0
			ref_idx = 0

			# access reference and alternate depths, and 
			# calculate allele balance
			rd = v.format('AD')[:,ref_idx]
			ad = v.format('AD')[:,alt_idx + 1]
			td = rd + ad
			ab = ad / td

			# access genotype qualities
			gq = v.gt_quals

			# get a list of sample indices that meet quality thresholds 
			# e.g., non-UNK, DP>=10, GQ>=20
			good_idxs = get_good_idxs(gts_reformatted, gq, td, 
										min_dp = args.min_dp, min_gq = args.min_gq)
			
			if good_idxs.shape[0] == 0: continue

			# make sure both founders have confident genotypes, so that we
			# can be sure that potential "singletons" are not inherited
			founder_idxs = np.array([smp2idx[B6], smp2idx[D2]])
			if np.intersect1d(founder_idxs, good_idxs).shape[0] != 2: continue

			# make sure both founders are HOM_REF
			if np.sum(gts_reformatted[founder_idxs]) != 0: continue
			
			n_mice_with_var = np.sum((gts_reformatted[good_idxs] == 2) | (gts_reformatted[good_idxs] == 1))

			if n_mice_with_var < 2: continue

			# get the mutation context for this sites using mutyper if it's a SNP
			# otherwise, we'll just report the kmer is an "indel"
			kmer = ancestor.mutation_type(v.CHROM, v.start, ref_allele, alt_allele)
			if None in kmer: kmer = "N>N"
			else: kmer = '>'.join(list(kmer))
			
			# annotate variant with phastCons score
			con_score = -1
			records = tb.query(v.CHROM, v.POS, v.POS + 1)
			for r in records:
				con_score = float(r[-1])

			outvals = []

			# output samples with variants
			for idx in good_idxs:
				mouse = idx2smp[idx]
				founder_hap = None
				hap_fh_pref = '/'.join(args.haplotypes[0].split('/')[:-1]) + '/'
				mouse_hap_fh = hap_fh_pref + '{}_haplotypes.csv'.format(mouse)
				if not os.path.isfile(mouse_hap_fh): continue
				hap_tree = make_interval_tree(mouse_hap_fh, datacol=True, delim=',')
				founder_hap = hap_tree[str(v.CHROM)].search(v.start, v.end)

				if founder_hap is None or len(founder_hap) == 0: founder_hap = 'U'
				else: founder_hap = str(founder_hap[0].data)

				m_genotype = gts_reformatted[idx]
				
				if m_genotype == 0: continue
				if ab[idx] < 0.9: continue

				outvals.append([v.CHROM, v.start, v.end, mouse, kmer, founder_hap, 
							gts_reformatted[idx], td[idx], ab[idx], con_score])
			if len(outvals) < 2: continue
				
			for ov in outvals:
				print (','.join(list(map(str, ov))), file=outcsvfh)
			
if __name__ == "__main__":
	doctest.testmod()
	p = argparse.ArgumentParser()
	p.add_argument('--vcf', required=True, 
						help="path to BXD VCF")
	p.add_argument('--ref', required=True, 
						help="path to reference genome")
	p.add_argument('--chrom', required=True, 
						help='chromosome to limit analysis')
	p.add_argument('--haplotypes', required=True, nargs="*",
						help='list of HMM-derived haplotypes for all samples')
	p.add_argument('--pcons', required=True, 
						help='path to phastCons file for the chromosome being queried')
	p.add_argument('--out', required=True, 
						help='name of output file')
	p.add_argument('-min_dp', default=10, type=int, 
						help='minimum depth required of singletons (default = 10)')
	p.add_argument('-min_gq', default=20, type=int, 
						help='minimum genotype quality required of singletons (default = 20)')
	p.add_argument('-nmer', default=3, type=int, 
						help='length of k-mers we want to catalog (default = 3)')
	p.add_argument('-exclude', 
						help='path to BED of regions to exclude')
	args = p.parse_args()
	run(args)

