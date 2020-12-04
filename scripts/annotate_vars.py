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

def combine_chr_df(path, header=None, skiprows=None):
    main_df = None
    for f in glob.glob(path):
        if main_df is None:
            if skiprows is not None:
                main_df = pd.read_csv(f, header=header, skiprows=[skiprows])
            else:
                main_df = pd.read_csv(f, header=header)

        else:
            if skiprows is not None:
                main_df = pd.concat([main_df, pd.read_csv(f, header=header, skiprows=[skiprows])])
            else:
                main_df = pd.concat([main_df, pd.read_csv(f, header=header)])
    return main_df

def format_bxd_name(strain, strain_summary="/Users/tomsasani/Downloads/strain_summary.xlsx"):
    """
    since most BXD strains are formatted like
    BXD002_TwJ_NNNNN, we need to hack a little bit
    to format the strain ID properly for QTL later.

    strain: str() of strain ID

    """
    
    strain_summary = pd.read_excel(strain_summary)

    official2expanded = dict(zip(strain_summary['Official name'], 
                                 strain_summary['Expanded name']))

    # extract the BXDNNN string, where
    # NNN is the RIL number
    if 'BXD' not in strain:
        return strain
    bxd_line = strain.split('_')[0]

    # one special exception
    if bxd_line == 'BXD24': return 'BXD024/TyJ-Cep290rd16/J'
   
    # extract the BXD RIL number (the NNN from before)
    # and strip any preceding zeros
    temp = re.compile("([a-zA-Z]+)([0-9]+)([a-zA-Z]*)")
    res = temp.match(bxd_line).groups()
    bxd_num = res[1].lstrip('0')
   
    # some strains have trailing information after
    # the RIL number, like RwwJ or TyJ
    try:
        addtl = strain.split('_')[1]
    except IndexError: addtl = ''

    suff = res[2]
    
    formatted_bxd = 'BXD' + str(bxd_num) + suff
    
    off_name = strain
    
    matches = []
    
    for official_n in strain_summary['Official name']:
        if formatted_bxd in official_n:
            matches.append(official_n)

    matched_off_name = None
    
    # if there's an obvious match, return it
    if len(matches) == 1: 
        matched_off_name = matches[0]
    # if there's more than one possible match...
    elif len(matches) > 1:
        matched = False
        for m in matches:
            # check and see if the additional sample name
            # info is enough for a match
            if addtl in m and (formatted_bxd == m.split('/')[0] or \
                               formatted_bxd == m.split('-')[0]): 
                matched = True
                matched_off_name = m
        if not matched:
            for m in matches:
                if formatted_bxd == m.split('/')[0]:
                    matched = True
                    matched_off_name = m
        if not matched: 
            for m in matches:
            # check and see if the additional sample name
            # info is enough for a match
                if addtl in m: 
                    matched = True
                    matched_off_name = m
        # if not, probably a situation where the sample
        # manifest has two entries and the VCF has only one
        if not matched: matched_off_name = matches[0]

    return official2expanded[matched_off_name]
    
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

    path = '/Users/tomsasani/harrislab/bxd/hmm_haplotypes/{}_{}_haplotypes.csv'.format(sample, chrom)


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
p.add_argument("--strain_metadata")
p.add_argument("--var_dir")
p.add_argument("--callable_bp")
p.add_argument("--out")
args = p.parse_args()

### ---
# read in strain metadata and construct relevant dictionaries that
# map sample names to metadata
### ---

summary = pd.read_excel(args.strain_metadata)

# convert verbose filial generations to integer numbers of generations in the file
summary['gen_at_seq'] = summary['Generation at sequencing'].apply(get_generation,
                                                                  remove_backcrossed=True)

# map expanded strain names to generation numbers
#strain2inbreed_gen = dict(zip(summary['Expanded name'], summary['gen_at_seq']))
strain2inbreed_gen = dict(zip(summary['bam_name'], summary['gen_at_seq']))

# map expanded strain names to epochs
#strain2epoch = dict(zip(summary['Expanded name'], summary['Epoch']))
strain2epoch = dict(zip(summary['bam_name'], summary['Epoch']))

# map expanded strain names to genenetwork names
#strain2genenet = dict(zip(summary['Expanded name'], summary['GeneNetwork name']))
strain2genenet = dict(zip(summary['bam_name'], summary['GeneNetwork name']))

# map expanded strain names to the number of generations
# each strain was intercrossed prior to inbreeding (0 for F2-derived strains)
strain2intercross_gens = dict(zip(summary['bam_name'],
                                  summary['n_advanced_intercross_gens']))

# map strain names to callable base pairs (referred to as "denominators")
denominators = pd.read_csv(args.callable_bp, names=['samp', 'denominator'])
strain2denom = dict(zip(denominators['samp'], denominators['denominator']))


### ---
# read in the singleton file, containing one line per singleton
### ---

singleton = combine_chr_df(args.var_dir + "*.exclude.csv", skiprows=0)

# rename columns in the combined dataframe
singleton.columns = ['chrom', 'start', 'end', 'bam_name', 'kmer', 'haplotype',
                     'gt', 'dp', 'ab', 'phastCons']

# remove very low and very high-depth mutations
singleton = singleton.query('dp >= 10 & dp < 100')

# add a column describing the 1-mer mutation type corresponding to all
# 3-mer "kmers" in the dataframe
singleton['base_mut'] = singleton['kmer'].apply(to_base_mut, cpg=True)

# remove any potential singletons identified in founder genomes. we're
# only interested in singletons observed in BXD RILs
founder_samps = ['4512-JFI-0333_C57BL_6J_two_lanes_phased_possorted_bam',
                 '4512-JFI-0334_DBA_2J_three_lanes_phased_possorted_bam']

singleton = singleton[~singleton['bam_name'].isin(founder_samps)]

# reformat BXD strain names from BAM notation
singleton['bxd_strain_conv'] = singleton['bam_name'].apply(lambda x: convert_bxd_name(x))

# map each of these reformatted BXD strain names to its "Expanded name" as
# reported in the Excel metadata file
#uniq_bxd_names = pd.unique(singleton['bxd_strain_conv'])
#expanded_strain_names = [format_bxd_name(n) for n in uniq_bxd_names]
#orig2expanded = dict(zip(uniq_bxd_names, expanded_strain_names))

# then, add a column to the dataframe with the "Expanded name" for each strain
#singleton['expanded_name'] = singleton['bxd_strain_conv'].apply(lambda s: orig2expanded[s])

# remove co-isogenic samples
iso_samps = ['4512-JFI-0348_BXD24_TyJ_Cep290_J_phased_possorted_bam',
             '4512-JFI-0382_BXD048a_RwwJ_phased_possorted_bam',
             '4512-JFI-0387_BXD65a_RwwJ_phased_possorted_bam'
             '4512-JFI-0388_BXD65b_RwwJ_phased_possorted_bam',
             '4512-JFI-0439_BXD73a_RwwJ_phased_possorted_bam'
             '4512-JFI-0440_BXD073b_RwwJ_phased_possorted_bam']

singleton = singleton[~singleton['bam_name'].isin(iso_samps)]

# use the "Expanded name" for each BXD strain to add columns of metadata to the dataframe
#singleton['gene_network_name'] = singleton['expanded_name'].apply(lambda s: strain2genenet[s])
#singleton['epoch'] = singleton['expanded_name'].apply(lambda s: strain2epoch[s])
#singleton['n_inbreeding_gens'] = singleton['expanded_name'].apply(lambda s: strain2inbreed_gen[s])
#singleton['n_intercross_gens'] = singleton['expanded_name'].apply(lambda s: strain2intercross_gens[s])
#singleton['n_callable_bp'] = singleton['bxd_strain'].apply(lambda s: strain2denom[s] if s in strain2denom else -1)

singleton['epoch'] = singleton['bam_name'].apply(lambda s: strain2epoch[s])
singleton['n_inbreeding_gens'] = singleton['bam_name'].apply(lambda s: strain2inbreed_gen[s])
singleton['n_intercross_gens'] = singleton['bam_name'].apply(lambda s: strain2intercross_gens[s])
singleton['n_callable_bp'] = singleton['bam_name'].apply(lambda s: strain2denom[s] if s in strain2denom else -1)

singleton['haplotype_at_qtl'] = singleton['bam_name'].apply(lambda s: find_haplotype(s))

singleton.to_csv(args.out, index=False)
