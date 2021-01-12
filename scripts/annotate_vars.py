import matplotlib.pyplot as plt
import pandas as pd
import csv
import scipy.stats as ss
from collections import defaultdict, Counter
from quicksect import IntervalTree
import numpy as np
import seaborn as sns
import glob
import sys
import math
import argparse
import re

def to_base_mut(k, cpg=False):
    if k == '0' or k == 0: return '0'
    if k == 'indel': return 'indel'
    nuc_1, nuc_2 = k.split('>')
    if not (len(nuc_1) == 3 and len(nuc_2) == 3): return 'indel'
    base = 'na'
    if len(nuc_1) % 2 == 1:
        mid_idx = int(len(nuc_1) / 2)
        base = nuc_1[mid_idx] + '>' + nuc_2[mid_idx]
        subseq_nuc = nuc_1[mid_idx + 1]
        if cpg:
            if nuc_1[mid_idx] == 'C' and subseq_nuc == 'G' and nuc_2[mid_idx] == 'T': base = 'CpG>TpG'
    else: base = k
    return base

def combine_chr_df(path):
    main_df = None
    for f in glob.glob(path):
        if main_df is None:
            main_df = pd.read_csv(f)
        else:
            main_df = pd.concat([main_df, pd.read_csv(f)])
    return main_df

def convert_bxd_name(name: str) -> str:
    """
    depending on the VCF file we're using we may
    have to convert BXD strain names for portability
    """
    bxd_line_name = '_'.join(name.split('_phased')[0].split('_')[1:])
    bxd_line_num = name.split('_phased')[0].split('_')[0].split('-')[-1]

    bxd_line_new = bxd_line_name + '_' + bxd_line_num

    return bxd_line_new

def get_generation(gen, remove_backcrossed=False):
    split = None
    try:
        split = re.split('(\d+)', gen)
    except TypeError: return 'NA'

    cur_gen = 0

    for i,e in enumerate(split):
        if 'F' in e:
            cur_gen += int(split[i + 1])
        elif 'N' in e:
            if remove_backcrossed: return 'NA'
            else:
                cur_gen /= (2 ** int(split[i + 1]))
        else: continue

    return int(cur_gen)

def find_haplotype(sample, location="chr4:115000000-118000000"):

    chrom = location.split(':')[0]
    start, end = location.split(':')[-1].split('-')

    path_pref = "/Users/tomsasani/harrislab/bxd/hmm_haplotypes"

    path = "{}/{}_{}_haplotypes.csv".format(path_pref, sample, chrom)

    tree = defaultdict(IntervalTree)
    added = defaultdict(int)
    f = gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'r')
    fh = csv.reader(f, delimiter=',')
    for i, line in enumerate(fh):
        tree[line[0]].add(int(line[1]), int(line[2]), other=line[-1])

    hap_at_loc = tree[chrom].search(int(start), int(end))

    hap = None

    for interval in hap_at_loc:
        if interval.start < int(start) and interval.end > int(end):
            hap = interval.data

    return hap

p = argparse.ArgumentParser()
p.add_argument("--strain_metadata", required=True, 
                    help="""metadata about strains in Excel format.""")
p.add_argument("--var_dir", required=True,
                    help="""directory containing per-chromosome BED files \
                            of variants""")
p.add_argument("--callable_bp", required=True,
                    help="""file containing a column with sample IDs and a \
                            column with the number of haploid callable base \
                            pairs in those samples.""")
p.add_argument("--out", required=True,
                    help="""name of the output extended BED file containing \
                            annotated variants""")
args = p.parse_args()

# read in strain metadata and construct relevant dictionaries that
# map sample names to metadata
summary = pd.read_excel(args.strain_metadata)

# convert verbose filial generations to integer numbers of generations in the file
summary['gen_at_seq'] = summary['Generation at sequencing'].apply(get_generation,
                                                                  remove_backcrossed=True)

# map expanded strain names to generation numbers
strain2inbreed_gen = dict(zip(summary['bam_name'], summary['gen_at_seq']))

# map expanded strain names to epochs
strain2epoch = dict(zip(summary['bam_name'], summary['Epoch']))

# map expanded strain names to genenetwork names
strain2genenet = dict(zip(summary['bam_name'], summary['GeneNetwork name']))

# map expanded strain names to the number of generations
# each strain was intercrossed prior to inbreeding (0 for F2-derived strains)
strain2intercross_gens = dict(zip(summary['bam_name'],
                                  summary['n_advanced_intercross_gens']))

# map strain names to callable base pairs (referred to as "denominators")
denominators = pd.read_csv(args.callable_bp, names=['samp', 'denominator'])
strain2denom = dict(zip(denominators['samp'], denominators['denominator']))

# read in the file of variants, containing one BED-format line per variant
variants = combine_chr_df(args.var_dir + "*.exclude.csv")

# remove very low and very high-depth mutations
variants = variants.query('dp >= 10 & dp < 100')

variants = variants.query('ab >= 0.9')

# add a column describing the 1-mer mutation type corresponding to all
# 3-mer "kmers" in the dataframe
variants['base_mut'] = variants['kmer'].apply(to_base_mut, cpg=True)

# remove any potential variantss identified in founder genomes. we're
# only interested in variantss observed in BXD RILs
founder_samps = ['4512-JFI-0333_C57BL_6J_two_lanes_phased_possorted_bam',
                 '4512-JFI-0334_DBA_2J_three_lanes_phased_possorted_bam']

variants = variants[~variants['bxd_strain'].isin(founder_samps)]

print (len(pd.unique(variants['bxd_strain'])))

# remove co-isogenic samples
iso_samps = ['4512-JFI-0348_BXD24_TyJ_Cep290_J_phased_possorted_bam',
             '4512-JFI-0347_BXD024_TyJ_phased_possorted_bam',
             '4512-JFI-0345_BXD029_Tlr4_J_phased_possorted_bam',
             '4512-JFI-0344_BXD29_Ty_phased_possorted_bam',
             '4512-JFI-0355_BXD152_phased_possorted_bam',
             '4512-JFI-0362_BXD155_phased_possorted_bam',
             '4512-JFI-0482_BXD087_RwwJ_phased_possorted_bam',
             '4512-JFI-0485_BXD194_redo_phased_possorted_bam',
             '4512-JFI-0382_BXD048a_RwwJ_phased_possorted_bam',
             '4512-JFI-0387_BXD65a_RwwJ_phased_possorted_bam'
             '4512-JFI-0388_BXD65b_RwwJ_phased_possorted_bam',
             '4512-JFI-0439_BXD73a_RwwJ_phased_possorted_bam'
             '4512-JFI-0440_BXD073b_RwwJ_phased_possorted_bam']

variants = variants[~variants['bxd_strain'].isin(iso_samps)]


# reformat BXD strain names from BAM notation
variants['bxd_strain_conv'] = variants['bxd_strain'].apply(lambda x: convert_bxd_name(x))

# add columns to the dataframe with relevant metadata
variants['epoch'] = variants['bxd_strain'].apply(lambda s: strain2epoch[s])
variants = variants.query('epoch != 6')
variants['n_inbreeding_gens'] = variants['bxd_strain'].apply(lambda s: strain2inbreed_gen[s])
variants = variants[variants['n_inbreeding_gens'] != "NA"]
variants = variants.query('n_inbreeding_gens >= 20')
variants['n_intercross_gens'] = variants['bxd_strain'].apply(lambda s: strain2intercross_gens[s])
variants['n_callable_bp'] = variants['bxd_strain_conv'].apply(lambda s: strain2denom[s] if s in strain2denom else "NA")
variants['haplotype_at_qtl'] = variants['bxd_strain'].apply(lambda s: find_haplotype(s))

variants.to_csv(args.out, index=False)
