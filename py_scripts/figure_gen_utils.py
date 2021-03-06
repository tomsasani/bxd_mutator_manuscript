import pandas as pd
import re
from collections import Counter, defaultdict
from quicksect import IntervalTree
import gzip
import csv
from typing import List, Any
import numpy as np


def to_base_mut(k: str, cpg: Any = False) -> str:
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


def combine_chr_df(path_list: List[str]) -> pd.DataFrame:
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
    """

    # split generation designations by "+" symbols, which
    # indicate transitions between original facilities and JAX
    split = None
    try:
        split = re.split("(\d+)", gen)
    except TypeError:
        return 'NA'

    cur_gen = 0

    for i, e in enumerate(split):
        # add each number of filial generations to the
        # cumulative sum of generations
        if 'F' in e:
            cur_gen += int(split[i + 1])
        # "N"s in JAX designations indicate backcrossing
        # generations. we don't consider strains that have
        # been backcrossed, so return "NA"
        elif 'N' in e:
            return -1
        else:
            continue

    return int(cur_gen)


def find_haplotype(genos: pd.DataFrame, sample: str) -> str:
    """
    figure out whether each strain has a B or D haplotype,
    or is heterozygous, at the genotype marker at the peak
    of the QTL on chromosome 4
    """

    genos_in_smp = genos[sample].values
    geno_freq = Counter(genos_in_smp)
    total = sum([i[1] for i in geno_freq.items()])
    most_freq_geno = "H"
    for g in ["B", "D"]:
        if geno_freq[g] > (total * 0.5): most_freq_geno = g[0]
        else: continue

    return most_freq_geno


def make_interval_tree(
    path: str,
    datacol: Any = False,
    delim: str = '\t',
) -> defaultdict(IntervalTree):
    """
    generate an interval tree from a simple BED
    file with format chr:start:end
    """

    tree = defaultdict(IntervalTree)
    added = defaultdict(int)
    f = gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'r')
    fh = csv.reader(f, delimiter=delim)
    for i, line in enumerate(fh):
        if datacol:
            tree[line[0]].add(int(line[1]), int(line[2]), other=line[3])
        else:
            tree[line[0]].add(int(line[1]), int(line[2]))

    return tree


def revcomp(seq: str) -> str:
    """
    reverse complement a nucleotide sequence.
    """

    rc_nuc = {
        'A': 'T',
        'C': 'G',
        'T': 'A',
        'G': 'C',
    }

    seq_rev = seq[::-1]
    seq_rev_comp = ''.join([rc_nuc[n] for n in list(seq_rev)])

    return seq_rev_comp


def calculate_years_per_gen(row: pd.Series) -> float:
    """
    calculate generation times by simply dividing time elapsed
    since 2017 by the number of generations of inbreeding
    """
    gens = row['gen_at_seq']
    if gens == -1 or gens == "NA":
        return -1.
    else:
        years_of_breeding = 2017 - int(row['Year breeding started'])
        years_per_gen = years_of_breeding / int(gens)
        return years_per_gen


def convert_cosmic_mutation(row: pd.Series) -> str:
    """
    convert cosmic mutation notation to 3-mer 
    mutation notation so that it matches the BXD data
    """

    context = row['Subtype']
    mutation = row['Type']

    changed_from = mutation.split('>')[0]
    changed_to = mutation.split('>')[1]

    e_kmer = context[0] + changed_to + context[-1]

    if changed_from == "T":
        context = revcomp(context)
        e_kmer = revcomp(e_kmer)

    return context + '>' + e_kmer


def convert_toyko_mutation(sequence: str):
    """
    convert TOY-KO mutations (reported as a string of 50
    upstream nucleotides plus the mutation plus a string of 50
    downstream nucleotides) to notation that matches the BXD data
    """
    mutation = sequence.split('[')[-1].split(']')[0]
    left_flank_1bp = sequence.split('/')[0].split('[')[0][-1]
    right_flank_1bp = sequence.split('/')[-1].split(']')[-1][0]

    anc, der = mutation.split('/')

    kmer_anc = left_flank_1bp + anc + right_flank_1bp
    kmer_der = left_flank_1bp + der + right_flank_1bp

    if mutation not in ["C/A", "G/T"]: return 'not_CA'

    # reverse complement if necessary
    rc = False
    if mutation[0] == "G":
        rc = True

    if rc: return "{}>{}".format(revcomp(kmer_anc), revcomp(kmer_der))
    else: return "{}>{}".format(kmer_anc, kmer_der)


def find_groups(a: List[Any]) -> List[Any]:
    """
    function to get the indexes of shared "groups"
    in an arbitrary list
    """
    groups = []
    cur_val, last_idx = 0, 0
    for idx, val in enumerate(a):
        if idx == 0:
            cur_val = val
            continue
        if val != cur_val or idx == len(a) - 1:
            if idx == len(a) - 1:
                groups.append((last_idx, idx + 1, cur_val))
            else:
                groups.append((last_idx, idx, cur_val))
            last_idx = idx
            cur_val = val
        else:
            cur_val = val
    return groups


def calc_new_gens(n_gens: int) -> int:
    """
    use the method from Uchimura et al. (2015) to calculate
    the number of generations in which an observed homozgyous
    singleton could have occurred in a given strain.
    """
    p_k_vals = defaultdict(float)

    for k in np.arange(1, n_gens + 1):
        if k == 1:
            p_k = 0
            p_k_vals[k] = p_k
        elif k == 2:
            p_k = 0.25
            p_k_vals[k] = p_k
        else:
            p_k_1 = p_k_vals[k - 1]
            p_k_2 = p_k_vals[k - 2]

            p_k = (p_k_1 / 2) + (p_k_2 / 4)
            p_k_vals[k] = p_k

    l_n = 0
    for k in np.arange(1, n_gens + 1):
        l_n += ((n_gens - k) * p_k_vals[k])
    return l_n


def clr(X: np.ndarray) -> np.ndarray:
    """
    perform a centered log-ratio transform
    """
    # the geometric mean acts as the center of the composition
    geom_mean = np.power(np.prod(X, axis=1), 1 / X.shape[1])
    return np.log(X / geom_mean[:, None])
