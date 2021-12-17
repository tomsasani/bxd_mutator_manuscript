import pytest
import pandas as pd
import numpy as np
from py_scripts.figure_gen_utils import (
    calc_new_gens,
    calculate_years_per_gen,
    convert_bxd_name,
    convert_cosmic_mutation,
    convert_toyko_mutation,
    find_groups,
    to_base_mut,
    get_generation,
    revcomp,
    clr,
)


def test_convert_bxd_name():
    assert convert_bxd_name("4512-JFI-0469_BXD100_RwwJ_phased_possorted_bam") == "BXD100_RwwJ_0469"
    assert convert_bxd_name("4512-JFI-0410_BXD013_TyJ_phased_possorted_bam") == "BXD013_TyJ_0410"


@pytest.mark.parametrize('input_mut,output_mut,cpg', [
    ("ACA>ATA", "C>T", True),
    ("CCG>CTG", "CpG>TpG", True),
    ("ACG>ATG", "C>T", False),
])
def test_to_base_mut(input_mut, output_mut, cpg):
    assert to_base_mut(input_mut, cpg=cpg) == output_mut


def test_get_generation():
    assert get_generation("F48N5F75+10pF12") == -1
    assert get_generation("F44+F14pF22") == 80


def test_revcomp():
    assert revcomp('ATTCAG') == 'CTGAAT'
    assert revcomp('T') == 'A'
    assert revcomp('CGA') == 'TCG'

def test_get_generation():
    assert get_generation("F48N5F75+10pF12") == -1
    assert get_generation("F44+F14pF22") == 80


@pytest.mark.parametrize(
    'row,exp_output',
    [(pd.Series({
        'gen_at_seq': 56,
        'Year breeding started': 1976
    }), 0.7321428),
     (pd.Series({
         'gen_at_seq': 12,
         'Year breeding started': 1999
     }), 1.5)])
def test_calculate_years_per_gen(row, exp_output):
    assert np.allclose(calculate_years_per_gen(row), exp_output)


@pytest.mark.parametrize('row,exp_output', [
    (pd.Series({
        'Subtype': 'ACA',
        'Type': 'C>G'
    }), 'ACA>AGA'),
    (pd.Series({
        'Subtype': 'TCG',
        'Type': 'C>G'
    }), 'TCG>TGG'),
])
def test_convert_cosmic_mutation(row, exp_output):
    assert convert_cosmic_mutation(row) == exp_output


def test_convert_toyko_mutation(toy_toyko_mutation):
    assert convert_toyko_mutation(toy_toyko_mutation) == "TCA>TAA"

@pytest.mark.parametrize('sequence', [([0, 0, 0, 0, 1, 1, 2, 3, 3, 3])])
def test_find_groups(sequence):
    assert find_groups(sequence) == [(0, 4, 0), (4, 6, 1), (6, 7, 2), (7, 10, 3)]


@pytest.mark.parametrize('gens,exp_output', [
    (24, 18.037886),
    (56, 50.),
    (90, 84.),
])
def test_calc_new_gens(gens, exp_output):
    assert np.allclose(calc_new_gens(gens), exp_output)


def test_clr(clr_array):
    assert np.allclose(
        clr(clr_array),
        np.array([[
            [0.28364631, -0.28364631],
            [-0.45460406, 0.45460406],
        ]]))
