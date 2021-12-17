import tabix
import numpy as np
import argparse
from collections import defaultdict
from cyvcf2 import VCF
from mutyper.ancestor import Ancestor
import doctest
import pandas as pd
from figure_gen_utils import make_interval_tree
from singleton_calling_utils import *
import time

def run(args):

    # -----------
    # read in VCF
    # -----------
    wild = "https://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/vcf/AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz"
    vcf = VCF(wild, gts012=True)

    # ------------
    # define some global variables that will come in handy
    # ------------
    UNK, HOM_REF, HET, HOM_ALT = range(-1, 3)

    ancestor = Ancestor(args.ref, k=args.nmer, sequence_always_upper=True)

    # get a list of all possible 3-mer mutation types
    mutations = enumerate_mutations(args.nmer)

    mut2idx = dict(zip(mutations, range(len(mutations))))

    samples = vcf.samples
    samples2use = [s for s in samples if s.startswith("Mmd")]

    vcf.set_samples(samples2use)

    # map sample indices to sample IDs and vice versa
    idx2smp = dict(zip(range(len(vcf.samples)), vcf.samples))
    smp2idx = dict(zip(vcf.samples, range(len(vcf.samples))))

    # generate an output array of size n_samples x n_mutations,
    # where we'll store mutation counts
    out_a = np.zeros((len(smp2idx), len(mut2idx)), dtype=np.int64)

    # map "ancestries" (i.e., Mmd subspecies) to sample indexes
    anc2idx = defaultdict(list)
    anc = [s.split('_')[1][:3] for s in vcf.samples]
    smp2anc = dict(zip(vcf.samples, anc))
    for smp in smp2idx:
        anc = smp2anc[smp]
        idx = smp2idx[smp]
        anc2idx[anc].append(idx)

    smp2founder_allele = np.zeros((2, len(smp2idx)), dtype=np.int64)

    # ------------
    # iterate over the VCF
    # ------------

    vcf_h = vcf(args.chrom if 'chr' in args.chrom else f'chr{args.chrom}')

    # generate an exclude file if it's provided
    exclude = None
    if args.exclude:
        exclude = make_interval_tree(args.exclude)

    # read in fixed mutations in D2 and B6
    fixed_vars = pd.read_csv(args.fixed_vars)

    start = time.time()
    for i, v in enumerate(vcf_h):

        if i % 50000 == 0 and i > 0:
            print(v.CHROM, v.POS, i / (time.time() - start))
            start = time.time()

        # check if the current variant is in the exclude file
        if exclude and len(exclude[v.CHROM].search(v.start, v.end)) > 0:
            continue
        if v.var_type != "snp": continue
        if '*' in v.ALT: continue
        if len(v.REF) > 1: continue
        if len(v.ALT) > 1: continue

        d2_allele, b6_allele = None, None

        if v.POS in fixed_vars['pos'].values:
            fixed_var_idx = np.where(fixed_vars['pos'].values == v.POS)
            d2_allele = fixed_vars['d2_allele'].values[fixed_var_idx]
            b6_allele = fixed_vars['b6_allele'].values[fixed_var_idx]

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
            ab = ad / td
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

            if d2_allele is not None and b6_allele is not None:
                b6_match = np.where(gts_reformatted == b6_allele)
                d2_match = np.where(gts_reformatted == d2_allele)
                smp2founder_allele[0][b6_match] += 1
                smp2founder_allele[1][d2_match] += 1

            smps_with_mut = np.where((gts_reformatted == 1)
                                     | (gts_reformatted == 2))[0]
            good_smps_with_mut = np.intersect1d(smps_with_mut, good_idxs)

            # skip this site if more than five samples
            # have a variant
            if np.sum(gts_reformatted > 0) > 5: continue

            out_a[good_smps_with_mut, kmer_idx] += 1

    out_df = []

    for smp in smp2idx:
        s_i = smp2idx[smp]
        d2_count = smp2founder_allele[1][s_i]
        b6_count = smp2founder_allele[0][s_i]
        ca_count = out_a[s_i, mut2idx["C>A"]]
        total_mut_count = np.sum(out_a[s_i])

        out_df.append({
            'sample': smp,
            'chrom': args.chrom,
            'ancestry': smp2anc[smp],
            'd2_count': d2_count,
            'b6_count': b6_count,
            'ca_count': ca_count,
            'mut_count': total_mut_count,
        })

    df = pd.DataFrame(out_df)

    df.to_csv(args.out)


if __name__ == "__main__":
    doctest.testmod()
    p = argparse.ArgumentParser()
    p.add_argument(
        '--ref',
        required=True,
        help="path to reference genome",
    )
    p.add_argument(
        "--fixed_vars",
        required=True,
        help='path to file with fixed differences between D and B',
    )
    p.add_argument("--chrom")
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
