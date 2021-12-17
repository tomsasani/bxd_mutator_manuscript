import numpy as np
from quicksect import IntervalTree
from typing import Tuple
import itertools
from figure_gen_utils import revcomp


def reformat_genotypes(
    gts: np.ndarray,
    alt_gt: int = 1,
) -> np.ndarray:
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


def enumerate_mutations(k):
    """
    generate a list of all possible kmer mutations
    """

    nmers = [''.join(x) for x in itertools.product('ATCG', repeat=k)]

    expanded_muts = []

    for n1 in nmers:
        for n2 in nmers:
            mut = '>'.join([n1, n2])
            if mut in expanded_muts: continue
            expanded_muts.append('>'.join([n1, n2]))

    expanded_rc = []
    for mut in expanded_muts:
        anc, der = mut.split('>')
        middle_nuc_anc = anc[int(len(anc) / 2)]
        middle_nuc_der = der[int(len(anc) / 2)]
        if middle_nuc_anc == middle_nuc_der: continue
        if levenshtein(anc, der) != 1: continue
        if middle_nuc_anc not in ('C', 'A'):
            anc, der = revcomp(anc), revcomp(der)
        expanded_rc.append('>'.join([anc, der]))

    expanded_rc_uniq = []
    for mut in expanded_rc:
        if mut in expanded_rc_uniq: continue
        expanded_rc_uniq.append(mut)

    return expanded_rc_uniq


def levenshtein(s1: str, s2: str):
    """
    calculate edit distance between two sequences
    """

    if len(s1) < len(s2):
        return levenshtein(s2, s1)

    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row

    return previous_row[-1]


def get_singleton_idx(
    gts: np.ndarray,
    ad: np.ndarray,
    rd: np.ndarray,
    ab_thresh: float = 0.75,
) -> int:
    """
	return the index of a sample with a putative
	singleton variant 

	gts: np.array() of sample genotypes
	ad: np.array() of sample alternate depths
	rd: np.array() of sample reference depths
	ab_thresh: lower bound on allowed allele balance at HETs
	"""

    UNK, HOM_REF, HET, HOM_ALT = range(-1, 3)

    # get all unique GTs and their frequencies at this site
    unique, counts = np.unique(gts, return_counts=True)
    uc = dict(zip(unique, counts))

    gts_to_use = (HOM_ALT, HET)

    # if there is more than one sample with a HET or HOM_ALT
    # genotype here, it's not a HQ singleton
    het_ha_sum = 0
    if HET in uc: het_ha_sum += uc[HET]
    if HOM_ALT in uc: het_ha_sum += uc[HOM_ALT]
    if het_ha_sum > 1: return None

    # if the singleton is HOM_ALT, return it. we'll apply
    # filters to individual singletons later in the main script
    if HOM_ALT in uc and HOM_ALT in gts_to_use:
        return np.where(gts == HOM_ALT)[0][0]
    # if the singleton is HET, we need to apply
    # filtering on allele balance.
    elif HET in uc and HET in gts_to_use:
        gt_idx = np.where(gts == HET)[0][0]
        ab = ad[gt_idx] / float(ad[gt_idx] + rd[gt_idx])
        if ab >= ab_thresh: return gt_idx
        else: return None

