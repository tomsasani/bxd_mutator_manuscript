import numpy as np
import argparse
from collections import defaultdict
from cyvcf2 import VCF
from mutyper.ancestor import Ancestor
import doctest
import pandas as pd
from singleton_calling_utils import (
    enumerate_mutations,
    make_interval_tree,
    normalize_var,
    reformat_genotypes,
    get_good_idxs,
)


def run(args):

    # -----------
    # read in VCF
    # -----------
    WILD_VCFFILE = "http://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/vcf/AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz"
    vcf = VCF(WILD_VCFFILE, gts012=True)

    # ------------
    # define some global variables that will come in handy
    # ------------
    UNK, HOM_REF, HET, HOM_ALT = range(-1, 3)

    ancestor = Ancestor(
        args.ref,
        k=args.nmer,
        sequence_always_upper=True,
    )

    # get a list of all possible 3-mer mutation types
    mutations = enumerate_mutations(args.nmer)

    mut2idx = dict(zip(mutations, range(len(mutations))))

    # map sample indices to sample IDs and vice versa
    idx2smp = dict(zip(range(len(vcf.samples)), vcf.samples))
    smp2idx = dict(zip(vcf.samples, range(len(vcf.samples))))

    # generate an output array of size n_samples x n_mutations,
    # where we'll store mutation counts
    out_a = np.zeros((len(smp2idx), len(mut2idx)), dtype=np.int64)

    # map "ancestries" (i.e., mouse subspecies) to sample indexes
    anc2idx = defaultdict(list)
    anc = [s.split('_')[0] for s in vcf.samples]
    smp2anc = dict(zip(vcf.samples, anc))
    for smp in smp2idx:
        anc = smp2anc[smp]
        idx = smp2idx[smp]
        anc2idx[anc].append(idx)

    # ------------
    # iterate over the VCF
    # ------------
    vcf_h = vcf
    if args.region:
        vcf_h = vcf(args.region)

    # generate an exclude file if it's provided
    exclude = None
    if args.exclude:
        exclude = make_interval_tree(args.exclude)

    region_s, region_e = args.region.split(':')[-1].split('-')
    region_s, region_e = int(region_s), int(region_e)

    window_size = 50000
    windows = np.arange(region_s, region_e, window_size)

    out_df = []

    for region_start in windows:
        # define numpy array to store sample x mutation matrix
        # in this window
        out_a = np.zeros((len(smp2idx), len(mut2idx)), dtype=np.int64)
        window = "chr4:{}-{}".format(region_start, region_start + window_size)
        vcf_h = vcf(window)
        for i, v in enumerate(vcf_h):
            # check if the current variant is in the exclude file
            if exclude and len(exclude[v.CHROM].search(v.start, v.end)) > 0:
                continue
            if v.var_type != "snp": continue
            if '*' in v.ALT: continue
            if len(v.REF) > 1: continue
            if len(v.ALT) > 1: continue
            # require all samples to have a callable genotype?
            #if v.call_rate < 1: continue

            for alt_idx, alt in enumerate(v.ALT):
                if len(alt) > 1: continue
                ref_allele, alt_allele = v.REF, alt
                # don't consider indels
                if len(v.REF) != len(alt): continue
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
                gts = np.array(v.genotypes)[:, :-1]
                gts_reformatted = reformat_genotypes(gts, alt_gt=alt_gt)

                # the index of the reference allele is always 0
                ref_idx = 0
                # access ref and alt depths
                rd = v.format('AD')[:, ref_idx]
                ad = v.format('AD')[:, alt_idx + 1]
                rd[rd < 0] = 0
                ad[ad < 0] = 0

                td = rd + ad
                gq = v.gt_quals

                # get a list of sample indices that meet quality thresholds
                # e.g., non-UNK, DP>=10, GQ>=20
                good_idxs = get_good_idxs(
                    gts_reformatted,
                    gq,
                    td,
                    min_dp=args.min_dp,
                    min_gq=args.min_gq,
                )

                if good_idxs.shape[0] == 0: continue

                # get the mutation context for this sites using mutyper if it's a SNP
                # otherwise, we'll just report the kmer is an "indel"
                kmer = ancestor.mutation_type(v.CHROM, v.start, ref_allele,
                                              alt_allele)
                if None in kmer: kmer = "N>N"
                else: kmer = '>'.join(list(kmer))
                if kmer == "N>N": continue
                kmer_idx = mut2idx[kmer]

                smps_with_mut = np.where((gts_reformatted == 1)
                                         | (gts_reformatted == 2))[0]
                good_smps_with_mut = np.intersect1d(smps_with_mut, good_idxs)

                # determine the number of samples in each subspecies with
                # a high-quality ALT genotype
                anc2gts = np.zeros(len(anc2idx))
                ii = 0
                for anc in anc2idx:
                    anc_idxs = np.array(anc2idx[anc])
                    anc_good_idxs = np.intersect1d(good_smps_with_mut,
                                                   anc_idxs)
                    anc_gts = gts_reformatted[anc_good_idxs]
                    anc2gts[ii] = np.sum(anc_gts)
                    ii += 1

                # skip this site if more than one sample in a
                # subspecies has an ALT genotype
                # i.e., only look at singletons
                if np.sum(anc2gts) == 0: continue
                #if np.sum(anc2gts) > 5: continue
                if not np.all(anc2gts <= 1): continue

                out_a[good_smps_with_mut, kmer_idx] += 1

        # loop over sample x mutation matrix and
        # output data for this window for each sample
        for s_i in np.arange(out_a.shape[0]):
            for kmer_idx in np.arange(out_a.shape[1]):
                samp_name = idx2smp[s_i]
                mut_type = mutations[kmer_idx]
                count = out_a[s_i, kmer_idx]
                total = np.sum(out_a[s_i])
                frac = count / total
                out_df.append({
                    "kmer": mut_type,
                    "region": window,
                    "count": count,
                    "frac": frac,
                    "sample": samp_name,
                    "total": total,
                })

    df = pd.DataFrame(out_df)

    df.to_csv(args.out)


if __name__ == "__main__":
    doctest.testmod()
    p = argparse.ArgumentParser()
    p.add_argument('--vcf')
    p.add_argument(
        '--ref',
        required=True,
        help="path to reference genome",
    )
    p.add_argument(
        '--region',
        required=True,
        help='region to limit analysis',
    )
    p.add_argument(
        '--out',
        required=True,
        help='name of output file',
    )
    p.add_argument(
        '-min_dp',
        default=10,
        type=int,
        help='minimum depth required of singletons (default = 10)',
    )
    p.add_argument(
        '-min_gq',
        default=20,
        type=int,
        help='minimum genotype quality required of singletons (default = 20)',
    )
    p.add_argument(
        '-nmer',
        default=3,
        type=int,
        help='length of k-mers we want to catalog (default = 3)',
    )
    p.add_argument(
        '-exclude',
        help='path to BED of regions to exclude',
    )
    args = p.parse_args()
    run(args)
