import numpy as np
from collections import defaultdict
from quicksect import IntervalTree
import csv
import gzip
from typing import Tuple, Any


def reformat_genotypes(gts: np.ndarray, alt_gt: int = 1) -> np.ndarray:
    """
	reformat np.array() of genotypes in VCF so
	that we can query multi-allelic variants
	"""

    # access genotypes on the paternal and maternal haplotypes
    p_gts, m_gts = gts[:, 0], gts[:, 1]

    p_unk = np.where(p_gts == -1)[0]
    m_unk = np.where(m_gts == -1)[0]

    assert np.sum(p_unk == m_unk) == p_unk.shape[0]

    unk_gts = p_unk

    # boolean arrays indicating whether samples have the alternate
    # allele in question on each haplotype. each sample gets a 0
    # if they don't have the ALT allele on that haplotype, and a 1
    # if they do.
    p_gts_with_var = np.array(p_gts == alt_gt, dtype=np.int8)
    m_gts_with_var = np.array(m_gts == alt_gt, dtype=np.int8)

    # reformat genotypes from the boolean arrays above to resemble the
    # genotypes we'd access using `v.gt_types`. that is, samples get a
    # 0 if they're HOMREF, 1 if they're HET, 2 if they're HOM_ALT.
    gts_reformatted = p_gts_with_var + m_gts_with_var

    # reformat UNK gts to be -1
    gts_reformatted[unk_gts] = -1

    return gts_reformatted


def normalize_var(ref: str, alt: str) -> Tuple[str, str]:
    """
	method for normalizing a variant by hand

	>>> normalize_var('TCC', 'CCC')
	('T', 'C')
	>>> normalize_var('CAGGGGG', 'CTGGGGG')
	('A', 'T')
	>>> normalize_var('GAT', 'GAT')
	('N', 'N')
	"""

    ref_a = np.array(list(ref))
    alt_a = np.array(list(alt))

    diff_nucs = ref_a != alt_a
    diff_idxs = np.where(diff_nucs)[0]

    if diff_idxs.shape[0] == 0: return ('N', 'N')
    elif diff_idxs.shape[0] > 1: return ('N', 'N')
    else:
        return (ref_a[diff_idxs[0]], alt_a[diff_idxs[0]])


def get_good_idxs(
        gts: np.ndarray,
        gq: np.ndarray,
        td: np.ndarray,
        min_dp: int = 10,
        min_gq: int = 20,
) -> np.ndarray:
    """
	get a list of indices corresponding to samples
	that meet a set of reasonable filters.

	gts: 	array of sample genotypes
	gq: 	array of sample genotype qualities 
	td: 	array of sample depths
	
	>>> get_good_idxs(np.array([0, 0, 2, 0]), np.array([20, 4, 99, 33]), np.array([15, 9, 22, 4]))
	array([0, 2])
	>>> get_good_idxs(np.array([0, -1, -1, 1]), np.array([20, -1, -1, 47]), np.array([15, -1, -1, 23]))
	array([0, 3])
	>>> get_good_idxs(np.array([0, 0, 0, 0]), np.array([55, 99, 10, 34]), np.array([5, 9, 23, 56]))
	array([3])
	"""

    UNK, HOM_REF, HET, HOM_ALT = range(-1, 3)

    # get indices of non-unk genotypes
    called_gts = np.where(gts != UNK)[0]
    # and indices of genotypes that meet GQ and DP filters
    good_gq = np.where(gq >= min_gq)
    good_depth = np.where(td >= min_dp)

    # intersect all of the above indices together to get
    # the output indices that meet all filters
    good_sites = np.intersect1d(good_gq, good_depth)
    good_sites_not_unk = np.intersect1d(good_sites, called_gts)

    return good_sites_not_unk


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


def convert_bxd_name(name: str) -> str:
    """
	depending on the VCF file we're using we may
	have to convert BXD strain names for portability
	"""
    bxd_line_name = '_'.join(name.split('_phased')[0].split('_')[1:])
    bxd_line_num = name.split('_phased')[0].split('_')[0].split('-')[-1]

    bxd_line_new = bxd_line_name + '_' + bxd_line_num

    return bxd_line_new
