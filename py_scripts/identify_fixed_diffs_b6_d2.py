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
    vcf = VCF("ftp://ftp-mouse.sanger.ac.uk/current_snps/mgp.v5.merged.snps_all.dbSNP142.vcf.gz", gts012=True)
    #vcf = VCF(args.vcf, gts012=True)

    # ------------
    # define some global variables that will come in handy
    # ------------
    HOM_REF, HET, HOM_ALT, UNK = range(4)

    B6 = 'C57BL_6NJ'
    D2 = 'DBA_2J'

    # -----------
    # generate a few necessary dictionaries
    # -----------
    
    vcf.set_samples([B6, D2])

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

    csvheader = ['chrom', 'pos', 'd2_allele', 'b6_allele']

    print (','.join(csvheader), file=outcsvfh)

    # ------------
    # iterate over the VCF
    # ------------
    # generate an exclude file if it's provided
    exclude = None
    if args.exclude:
        exclude = make_interval_tree(args.exclude)
    
    for i,v in enumerate(vcf(args.chrom)):
        if i % 50000 == 0 and i > 0: print (v.CHROM, v.POS)
        # check if the current variant is in the exclude file
        if exclude and len(exclude[v.CHROM].search(v.start, v.end)) > 0: continue
        if '*' in v.ALT: continue
        # limit to SNPs
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

            gts_reformatted = reformat_genotypes(gts, alt_gt=alt_gt)

            # the index of the reference allele is always 0
            ref_idx = 0

            # access reference and alternate depths, and 
            # calculate allele balance
            # rd = v.format('AD')[:,ref_idx]
            # ad = v.format('AD')[:,alt_idx + 1]
            # td = rd + ad
            # ab = ad / td

            # # access genotype qualities
            # gq = v.gt_quals

            # # get a list of sample indices that meet quality thresholds 
            # # e.g., non-UNK, DP>=10, GQ>=20
            # good_idxs = get_good_idxs(gts_reformatted, gq, td, 
            #                             min_dp = args.min_dp, min_gq = args.min_gq)

            # if good_idxs.shape[0] == 0: continue

            # make sure both founders have confident genotypes
            #founder_idxs = np.array([smp2idx[B6], smp2idx[D2]])
            #if np.intersect1d(founder_idxs, good_idxs).shape[0] != 2: continue

            d2_idx, b6_idx = smp2idx[D2], smp2idx[B6]
            
            # make sure D2 and B6 are opposite homozygotes
            if not np.abs(gts_reformatted[d2_idx] - gts_reformatted[b6_idx]) == 2: continue

            outvals = [v.CHROM, v.POS, gts_reformatted[b6_idx], gts_reformatted[d2_idx]]

            print (','.join(list(map(str, outvals))), file=outcsvfh)

if __name__ == "__main__":
    doctest.testmod()
    p = argparse.ArgumentParser()
    p.add_argument('--ref', required=True, 
                        help="path to reference genome")
    p.add_argument('--chrom', required=True, 
                        help='chromosome to limit analysis')
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

