import pytest

from py_scripts.singleton_calling_utils import get_good_idxs, get_singleton_idx, levenshtein, reformat_genotypes, normalize_var
import numpy as np

def test_reformat_genotypes(genotype_array):
    assert np.array_equal(reformat_genotypes(genotype_array), np.array([0, 1, 1, 2]))

def test_reformat_genotypes(genotype_array_with_unk):
    assert np.array_equal(reformat_genotypes(genotype_array_with_unk), np.array([-1, 1, -1, 2]))


@pytest.mark.parametrize('invara,invarb,outvar', [
    ('TCC', 'CCC', ('T', 'C')),
    ('CAGGGGG', 'CTGGGGG', ('A', 'T')),
    ('GAT', 'GAT', ('N', 'N')),
])
def test_normalize_var(invara, invarb, outvar):
    assert normalize_var(invara, invarb) == outvar


@pytest.mark.parametrize('gt,dp,gq,idxs', [
    (
        np.array([0, 0, 2, 0]),
        np.array([20, 4, 99, 33]),
        np.array([15, 9, 22, 4]),
        np.array([0, 2]),
    ),
    (
        np.array([0, -1, -1, 1]),
        np.array([20, -1, -1, 47]),
        np.array([15, -1, -1, 23]),
        np.array([0, 3]),
    ),
    (
        np.array([0, 0, 0, 0]),
        np.array([55, 99, 10, 34]),
        np.array([5, 9, 23, 56]),
        np.array([3]),
    ),
])
def test_get_good_idxs(gt, dp, gq, idxs):
    assert np.array_equal(get_good_idxs(gt, dp, gq), idxs)


@pytest.mark.parametrize('seqa,seqb,dist', [
    ('AAA', 'ATA', 1),
    ('ATTTT', 'CAAAA', 5),
])
def test_levenshtein(seqa, seqb, dist):
    assert levenshtein(seqa, seqb) == dist


@pytest.mark.parametrize('gt,rd,ad,idx', [
    (
        np.array([0, 0, 2, 0]),
        np.array([0, 0, 12, 0]),
        np.array([15, 9, 0, 18]),
        2,
    ),
    (
        np.array([0, 0, 2, 1]),
        np.array([0, 0, 12, 4]),
        np.array([15, 9, 0, 8]),
        None,
    ),
    (
        np.array([0, 0, 0, 1]),
        np.array([0, 0, 6, 10]),
        np.array([15, 9, 1, 0]),
        3,
    ),
])
def test_get_singleton_idx(gt, rd, ad, idx):
    assert get_singleton_idx(gt, rd, ad) == idx