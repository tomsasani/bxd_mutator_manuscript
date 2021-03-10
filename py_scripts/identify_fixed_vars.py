import tabix
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

	# name the output CSV file, which will also contain all
	# passing sites
	outcsvfh = open(args.out, 'w')

	csvheader = ['chrom', 'start', 'end', 'bxd_strain', 'kmer', 
				'haplotype', 'gt', 'dp', 'ab', 'phastCons']

	print (','.join(csvheader), file=outcsvfh)

	# ------------
	# iterate over the VCF
	# ------------
	# generate an exclude file if it's provided
	exclude = None
	if args.exclude:
		exclude = make_interval_tree(args.exclude)

	# loop over each 50kbp window across the genome
	with open(args.intervals, 'r') as interval_fh:
		intervals = csv.reader(interval_fh, delimiter='\t')
		for interval in intervals:
			i_chrom, i_start, i_end = interval
			# limit to intervals on the chromosome of interest
			if i_chrom != args.chrom: continue
			# to avoid counting every mutation in the VCF,
			# we'll only count N mutations per interval
			mutations_counted = 0
			ival = '{}:{}-{}'.format(i_chrom, i_start, i_end)
			for i,v in enumerate(vcf(ival)):
				# check if the current variant is in the exclude file
				if exclude and len(exclude[v.CHROM].search(v.start, v.end)) > 0: continue
				if '*' in v.ALT: continue
				# limit to SNPs
				if v.var_type != "snp": continue
				if len(v.REF) > 1: continue

				if mutations_counted > 5: break

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

					# make sure both founders have confident genotypes
					founder_idxs = np.array([smp2idx[B6], smp2idx[D2]])
					if np.intersect1d(founder_idxs, good_idxs).shape[0] != 2: continue
					
					# get the mutation context for this sites using mutyper if it's a SNP
					# otherwise, we'll just report the kmer is an "indel"
					kmer = ancestor.mutation_type(v.CHROM, v.start, ref_allele, alt_allele)
					if None in kmer: kmer = "N>N"
					else: kmer = '>'.join(list(kmer))
					
					# annotate variant with phastCons score
					con_score = -1
					tb = tabix.open(args.pcons)
					records = tb.query(v.CHROM, v.start, v.end)
					for r in records:
						con_score = float(r[-1])

					# output a single sample's genotype at this site
					s_idx = np.random.choice(good_idxs, 1)[0]
					# make sure sample has a "fixed" genotype
					if gts_reformatted[s_idx] != HOM_ALT: continue
					bxd_line = idx2smp[s_idx]

					# since these variants are fixed in D2, haplotype on which
					# the variant occurred is always D2 -- we don't actually
					# care which haplotype these occurred on, but the column needs
					# to be present for downstream processing
					haplotype = "D"
					
					outvals = [v.CHROM, v.start, v.end, bxd_line, kmer, haplotype,
							gts_reformatted[s_idx], td[s_idx], ab[s_idx], con_score]

					print (','.join(list(map(str, outvals))), file=outcsvfh)

					mutations_counted += 1

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
						help='list of paths to HMM-derived haplotypes for all samples')
	p.add_argument('--pcons', required=True, 
						help='path to phastCons file for the chromosome being queried')
	p.add_argument('--intervals', required=True, 
						help='file containing 50kbp windows in mm10')
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

