import tabix
import numpy as np
import gzip
import csv
import argparse
from cyvcf2 import VCF
from mutyper.ancestor import Ancestor
import doctest
from singleton_calling_utils import *

def get_singleton_idx(gts: np.array(int), 
					  ad: np.array(int), 
					  rd: np.array(int), 
					  ab_thresh=0.75) -> int:
	"""
	return the index of a sample with a putative
	singleton variant 

	gts: np.array() of sample genotypes
	ad: np.array() of sample alternate depths
	rd: np.array() of sample reference depths
	ab_thresh: lower bound on allowed allele balance at HETs

	>>> get_singleton_idx(np.array([0, 0, 2, 0]), np.array([0, 0, 12, 0]), np.array([15, 9, 0, 18]))
	2
	>>> get_singleton_idx(np.array([0, 0, 2, 1]), np.array([0, 0, 12, 4]), np.array([15, 9, 0, 8])) is None
	True
	>>> get_singleton_idx(np.array([0, 0, 0, 1]), np.array([0, 0, 6, 10]), np.array([15, 9, 1, 0]))
	3
	"""

	UNK, HOM_REF, HET, HOM_ALT = range(-1, 3)

	# get all unique GTs and their frequencies at this site
	unique, counts = np.unique(gts, return_counts=True)
	uc = dict(zip(unique, counts))

	gts_to_use = (HOM_ALT, HET)

	# if there is more than one sample with a HET or HOM_ALT
	# genotype here, it's not a HQ singleton
	het_ha_sum = 0
	if HET in uc: het_ha_sum += uc[HET]
	if HOM_ALT in uc: het_ha_sum += uc[HOM_ALT]
	if het_ha_sum > 1: return None

	# if the singleton is HOM_ALT, return it. we'll apply
	# filters to individual singletons later in the main script
	if HOM_ALT in uc and HOM_ALT in gts_to_use: 
		return np.where(gts == HOM_ALT)[0][0]
	# if the singleton is HET, we need to apply
	# filtering on allele balance. 
	elif HET in uc and HET in gts_to_use:
		gt_idx = np.where(gts == HET)[0][0]
		ab = ad[gt_idx] / float(ad[gt_idx] + rd[gt_idx])
		if ab >= ab_thresh: return gt_idx
		else: return None


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

			# get the mutation context for this sites using mutyper if it's a SNP
			# otherwise, we'll just report the kmer is an "indel"
			kmer = ancestor.mutation_type(v.CHROM, v.start, ref_allele, alt_allele)
			if None in kmer: kmer = "N>N"
			else: kmer = '>'.join(list(kmer))
			
			# find the singleton, allowing for high AB HETs
			singleton = get_singleton_idx(gts_reformatted, ad, rd, ab_thresh=0.9)
			if singleton is None: continue

			# make sure singleton also passes various filters
			if gq[singleton] < args.min_gq: continue
			if td[singleton] < args.min_dp: continue

			# get name of BXD with singleton
			bxd_line = idx2smp[singleton]
			if bxd_line in (B6, D2): continue

			# for the mouse with a derived singleton, figure out
			# which haplotype the variant is on by accessing the corresponding
			# interval tree from the BXD haplotype files
			singleton_founder_hap = None
			for haplotype_fh in args.haplotypes:
				if bxd_line not in haplotype_fh: continue
				hap_tree = make_interval_tree(haplotype_fh, datacol=True, delim=',')
				singleton_founder_hap = hap_tree[str(v.CHROM)].search(v.start, v.end)
			
			# make sure we've got a legit haplotype here
			if singleton_founder_hap is None or len(singleton_founder_hap) == 0: continue
			singleton_founder_hap = str(singleton_founder_hap[0].data)
			hap2hap = {"0": "B", "1": "D"}
			singleton_founder_hap = hap2hap[singleton_founder_hap]

			# make sure that there's at least one other mouse
			# with the same founder haplotype, but not the ALT
			samp_with_hap = 0
			samp_with_hap_hr = 0
			# loop over all BXD mice
			for mouse in samples_in_vcf:
				# skip the one with the putative singleont
				if mouse == bxd_line: continue
				# skip founders
				if mouse in (D2, B6): continue 
				# skip low-coverage founders
				if 'BXD' not in mouse: continue 
				founder_hap = None
				for haplotype_fh in args.haplotypes:
					if mouse not in haplotype_fh: continue
					hap_tree = make_interval_tree(haplotype_fh, datacol=True, delim=',')
					founder_hap = hap_tree[str(v.CHROM)].search(v.start, v.end)
				if founder_hap is None or len(founder_hap) == 0: continue
				founder_hap = str(founder_hap[0].data)
				founder_hap = hap2hap[founder_hap]
				if founder_hap != singleton_founder_hap: continue
				samp_with_hap += 1
				m_idx = smp2idx[mouse]
				if not gts_reformatted[m_idx] == HOM_REF: continue
				if not (gq[m_idx] >= args.min_gq and td[m_idx] >= args.min_dp): continue
				samp_with_hap_hr += 1

			# if we don't have another mouse with the same haplotype (but HOM_REF)
			# skip this site
			if samp_with_hap_hr < 1: continue

			# annotate variant with phastCons score
			con_score = -1
			records = tb.query(v.CHROM, v.start, v.end)
			for r in records:
				con_score = float(r[-1])

			outvals = [v.CHROM, v.start, v.end, bxd_line, kmer, 
						singleton_founder_hap, gts_reformatted[singleton], 
						td[singleton], ab[singleton], con_score]

			print (','.join(list(map(str, outvals))), file=outcsvfh)

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

