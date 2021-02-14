import doctest
import pandas as pd
import re
import glob
from collections import Counter

def to_base_mut(k: str, cpg=False) -> str:
    """
    convert a 3-mer mutation type into its 
    single-nucleotide mutation type

    >>> to_base_mut("ACA>ATA")
    "C>T"
    >>> to_base_mut("CCG>CTG", cpg=True)
    "CpG>TpG"
    >>> to_base_mut("ACG>ATG", cpg=False)
    "C>T"
    """

    if k == 'indel': return 'indel'

    # make sure the mutation type is a true 3-mer
    # (a 3-nucleotide reference to a 3-nucleotide alternate)
    nuc_1, nuc_2 = k.split('>')
    if not (len(nuc_1) == 3 and len(nuc_2) == 3): return 'indel'

    base = None

    # grab the middle nucleotide from either the reference
    # or alternate alleles
    mid_idx = 1
    base = nuc_1[mid_idx] + '>' + nuc_2[mid_idx]

    # special situation where mutation is CpG>TpG
    subseq_nuc = nuc_1[mid_idx + 1]
    if cpg:
        if base == "C>T" and subseq_nuc == 'G': 
            base = 'CpG>TpG'
        
    return base

def convert_bxd_name(name: str) -> str:
    """
    the names of BXDs are not consistent across all metadata 
    and VCF files, so we sometimes have to convert their names

    >>> convert_bxd_name("4512-JFI-0469_BXD100_RwwJ_phased_possorted_bam")
    "BXD100_RwwJ_0469"
    >>> convert_bxd_name("4512-JFI-0410_BXD013_TyJ_phased_possorted_bam")
    "BXD013_TyJ_0410"
    """
    bxd_line_name = '_'.join(name.split('_phased')[0].split('_')[1:])
    bxd_line_num = name.split('_phased')[0].split('_')[0].split('-')[-1]

    bxd_line_new = bxd_line_name + '_' + bxd_line_num

    return bxd_line_new

def combine_chr_df(path_list: []) -> pd.DataFrame():
    """
    singleton variants are stored in per-chromosome CSV files.
    before annotating or processing these variants, we first
    use this function to combine the per-chromosome CSVs into a
    single pandas DataFrame
    """
    main_df = None
    for f in path_list:
        if main_df is None:
            main_df = pd.read_csv(f)
        else:
            main_df = pd.concat([main_df, pd.read_csv(f)])
    return main_df

def get_generation(gen: str) -> int:
    """
    given a Jackson Labs-formatted string that designates
    the inbreeding/backcrossing history of a strain, calculate
    the total number of generations of inbreeding a strain has
    undergone, and make a note of strains that have been backcrossed

    >>> get_generation("F48N5F75+10pF12")
    -1
    >>> get_generation("F44+F14pF22")
    80
    """

    # split generation designations by "+" symbols, which
    # indicate transitions between original facilities and JAX
    split = None
    try:
        split = re.split('(\d+)', gen)
    except TypeError: 
        return 'NA'

    cur_gen = 0

    for i,e in enumerate(split):
        # add each number of filial generations to the 
        # cumulative sum of generations
        if 'F' in e:
            cur_gen += int(split[i + 1])
        # "N"s in JAX designations indicate backcrossing
        # generations. we don't consider strains that have
        # been backcrossed, so return "NA"
        elif 'N' in e: return -1
        else: continue

    return int(cur_gen)

def find_haplotype(genos: list, sample: str, rsids: []) -> str:
    """
    figure out whether each strain has a B or D haplotype,
    or is heterozygous, at the genotype marker at the peak
    of the QTL on chromosome 4
    """
    
    genos_in_smp = genos[sample].values
    geno_freq = Counter(genos_in_smp)
    most_freq_geno = "H"
    for g in ["B", "D"]:
        if geno_freq[g] > (len(rsids) * 0.75): most_freq_geno = g[0]
        else: continue
    
    return most_freq_geno

def revcomp(seq):
    """
    reverse complement a nucleotide sequence.

    >>> revcomp('ATTCAG')
    'CTGAAT'
    >>> revcomp('T')
    'A'
    >>> revcomp('CGA')
    'TCG'
    """

    rc_nuc = {'A':'T', 'C':'G', 'T':'A', 'G':'C'}

    seq_rev = seq[::-1]
    seq_rev_comp = ''.join([rc_nuc[n] for n in list(seq_rev)])

    return seq_rev_comp