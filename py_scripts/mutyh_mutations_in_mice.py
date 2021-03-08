import numpy as np
import argparse
from cyvcf2 import VCF
from singleton_calling_utils import get_good_idxs

def run(args):

	# -----------
	# read in VCF
	# -----------
	
	# htslib can read remote VCFs
	WILD_VCFFILE = "http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/vcf/AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz"
	MGP_VCFFILE = "ftp://ftp-mouse.sanger.ac.uk/current_snps/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"

	vcf = VCF(MGP_VCFFILE, gts012=True)
	if args.is_wild:
		vcf = VCF(WILD_VCFFILE, gts012=True)

	# ------------
	# define some global variables that will come in handy
	# ------------
	HOM_REF, HET, HOM_ALT = range(3)
	
	# -----------
	# generate a few necessary dictionaries
	# -----------
	# map sample indices to sample IDs and vice versa
	idx2smp = dict(zip(range(len(vcf.samples)), vcf.samples))
	smp2idx = dict(zip(vcf.samples, range(len(vcf.samples))))


	# define intervals corresponding to each of 6 Mutyh
	# missense mutations
	mutyh_vars = ["4:116814337-116814338",
				  "4:116814394-116814394",
				  "4:116815657-116815658",
				  "4:116817415-116817416",
				  "4:116817418-116817419",
				  "4:116816475-116816476"]
	if args.is_wild:
		mutyh_vars = ['chr' + i for i in mutyh_vars]

	outfh = open(args.out, 'w')

	# loop over intervals and catalog samples with mutations
	for interval in mutyh_vars:
		for i,v in enumerate(vcf(interval)):
			if v.FILTER not in (None, "PASS"): continue
			if v.var_type != "snp": continue
			
			# get arrays of sample genotypes
			gts = v.gt_types

			# if we're using the wild mouse VCF, we can access 
			# genotype qualities and depth information. otherwise
			# we'll just artificially set the GQ and DP info to be
			# at "passing" thresholds
			gq = np.full(gts.shape[0], 99)
			td = np.full(gts.shape[0], 99)
			
			if args.is_wild:
				gq = v.gt_quals
				ad, rd = v.gt_alt_depths, v.gt_ref_depths
				td = ad + rd

			# get a list of sample indices that meet quality thresholds 
			good_idxs = get_good_idxs(gts, gq, td)

			if good_idxs.shape[0] == 0: continue

			# loop over samples and print out their genotypes
			# at each of the mutyh missense mutations
			for smp in smp2idx:
				s_i = smp2idx[smp]
				genotype = -1
				if s_i in good_idxs:
					genotype = gts[s_i]
				print (','.join([interval, smp, str(genotype), 
									v.REF, v.ALT[0]]), file=outfh)

if __name__ == "__main__":
	p = argparse.ArgumentParser()
	p.add_argument('--out', required=True,
						help = 'path to outfile')
	p.add_argument('-is_wild', action='store_true',
						help = 'boolean flag to indicate whether the VCF is from wild mice')
	args = p.parse_args()
	run(args)

