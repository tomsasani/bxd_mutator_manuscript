import numpy as np
import argparse
import re
from collections import defaultdict, Counter
from cyvcf2 import VCF
from mutyper.ancestor import Ancestor
import itertools
import doctest
import pandas as pd
from scipy.stats import contingency
from singleton_calling_utils import *
import matplotlib.pyplot as plt
import seaborn as sns
from codon_table import codons
from pyfaidx import Fasta
import scipy.stats as ss

ref = Fasta("/Users/tomsasani/harrislab/bxd_mutator_ms/data/ref/mm10.fa")

codons = {k.replace('U', 'T'):v for k,v in codons.items()}

mutyh_exons = pd.read_csv("/Users/tomsasani/Downloads/mutyh.exons.gencode.txt", sep='\t').head(1)

exon_starts, exon_ends = mutyh_exons.exonStarts.values[0].split(','), mutyh_exons.exonEnds.values[0].split(',')

exon_positions = zip(exon_starts, exon_ends)

# -----------
# read in VCF
# -----------
vcf = VCF("/Users/tomsasani/harrislab/bxd_mutator_ms/vcf/wild.100.120.snpeff.vcf.gz", gts012=True)

# ------------
# define some global variables that will come in handy
# ------------
UNK, HOM_REF, HET, HOM_ALT = range(-1, 3)

# map sample indices to sample IDs and vice versa
idx2smp = dict(zip(range(len(vcf.samples)), vcf.samples))
smp2idx = dict(zip(vcf.samples, range(len(vcf.samples))))

# map "ancestries" (i.e., mouse subspecies) to sample indexes
anc2idx = defaultdict(list)
anc = [s.split('_')[0] for s in vcf.samples]
smp2anc = dict(zip(vcf.samples, anc))
for smp in smp2idx:
    anc = smp2anc[smp]
    idx = smp2idx[smp]
    anc2idx[anc].append(idx)

pos2mut = defaultdict()
pos2cons = defaultdict()

for s, e in list(exon_positions)[:-1]:
    
    window = "chr4:{}-{}".format(s, e)

    vcf_h = vcf(window)
    for v in vcf_h:
        if v.var_type != "snp": continue
        if '*' in v.ALT: continue
        if len(v.REF) > 1: continue
        if len(v.ALT) > 1: continue
        # require all samples to have a callable genotype?
        if v.call_rate < 1: continue

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
            gts = np.array(v.genotypes)[:,:-1]
            gts_reformatted = reformat_genotypes(gts, alt_gt=alt_gt)

            # the index of the reference allele is always 0
            ref_idx = 0
            # access ref and alt depths
            rd = v.format('AD')[:,ref_idx]
            ad = v.format('AD')[:,alt_idx + 1]
            rd[rd < 0] = 0
            ad[ad < 0] = 0

            td = rd + ad
            gq = v.gt_quals

            ab = ad / td

            high_ab = np.where(ab > 0.9)[0]
            gts_reformatted[high_ab] = 2

            # get a list of sample indices that meet quality thresholds
            # e.g., non-UNK, DP>=10, GQ>=20
            good_idxs = get_good_idxs(
                gts_reformatted,
                gq,
                td,
                min_dp=10,
                min_gq=1,
            )

            smps_with_mut = np.where((gts_reformatted == 1) | (gts_reformatted == 2))[0]
            good_smps_with_mut = np.intersect1d(smps_with_mut, good_idxs)

            if good_smps_with_mut.shape[0] == 0: continue
            #if good_idxs.shape[0] != len(smp2idx): continue

            # figure out whether any of the ancestries are fixed or polymorphic
            anc2state = defaultdict(str)
            ii = 0
            for anc in anc2idx:
                anc_idxs = np.array(anc2idx[anc])
                anc_good_idxs = np.intersect1d(good_idxs, anc_idxs)
                anc_gts = gts_reformatted[anc_good_idxs]
                if np.sum(anc_gts) == 0: 
                    anc2state[anc] = "no_mutation"
                elif np.sum(anc_gts) == anc_idxs.shape[0] * 2:
                    anc2state[anc] = "fixed"
                else:
                    anc2state[anc] = "polymorphic"
            
            if any([v == "polymorphic" for k,v in anc2state.items()]): 
                pos2mut[v.POS] = "polymorphic"
            elif all([(v == "fixed" or v == "no_mutation") for k,v in anc2state.items()]):
                pos2mut[v.POS] = "fixed"

            relevant_transcript = [ann for ann in v.INFO.get("ANN").split(',') if "ENSMUST00000102699" in ann]

            consequence = relevant_transcript[0].split('|')[1]
            pos2cons[v.POS] = consequence


fixed_s, fixed_ns = 0, 0
poly_s, poly_ns = 0, 0

for pos in pos2mut:
    if pos2mut[pos] == "polymorphic" and pos2cons[pos] == "synonymous_variant": poly_s += 1
    elif pos2mut[pos] == "polymorphic" and pos2cons[pos] == "missense_variant": poly_ns += 1
    elif pos2mut[pos] == "fixed" and pos2cons[pos] == "synonymous_variant": fixed_s += 1
    elif pos2mut[pos] == "fixed" and pos2cons[pos] == "missense_variant": fixed_ns += 1

contingency_table = [[fixed_s, fixed_ns], [poly_s, poly_ns]]
print (contingency_table)
print (ss.fisher_exact(contingency_table))


# https://nov2020.archive.ensembl.org/Mus_musculus/Transcript/Exons?db=core;g=ENSMUSG00000028687;r=4:116807723-116819431;t=ENSMUST00000102699

# https://nov2020.archive.ensembl.org/Mus_musculus/Transcript/Sequence_cDNA?db=core;g=ENSMUSG00000028687;r=4:116807723-116819431;t=ENSMUST00000102699
