import numpy as np
import pytest

@pytest.fixture
def toy_toyko_mutation():
    rng = np.random.default_rng(42)
    nucs = ['A', 'T', 'C', 'G']

    upstream_idxs = rng.integers(low=0, high=4, size=50, dtype=np.int8)
    downstream_idxs = rng.integers(low=0, high=4, size=50, dtype=np.int8)

    upstream_nucs = [nucs[i] for i in upstream_idxs]
    downstream_nucs = [nucs[i] for i in downstream_idxs]

    return ''.join(upstream_nucs) + f"[C/A]" + ''.join(downstream_nucs)


@pytest.fixture
def clr_array():
    rng = np.random.default_rng(42)
    orig = rng.random((1, 2))
    append = np.ones((1, 2)) - orig
    return np.vstack((orig, append))


@pytest.fixture
def genotype_array():
    return np.array([
        [0, 0, -1],
        [0, 1, -1],
        [0, 1, -1],
        [1, 1, -1],
    ])

@pytest.fixture
def genotype_array_with_unk():
    return np.array([
        [-1, -1, -1],
        [0, 1, -1],
        [-1, -1, -1],
        [1, 1, -1],
    ])
