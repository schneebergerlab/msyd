#!/bin/python3

import pansyri
import pansyri.classes.cigar as cigar
from pansyri.classes.cigar import Cigar
import pytest

@pytest.fixture(params=["10=", "10=2X4=", "100=20I4X5D10I40=", "10=30D5=2X2=20I10="])
def example_cigar():
    return Cigar.from_string("10=")

def test_from_strings():
    tuplelists = [ [ [8, '='], [2, 'I'], [4, '='], [2, 'D']], [[100, '=']], [[100, '='], [2, 'X'], [20, '='], [4, 'I'], [40, '='], [1, 'X']] ]
    strings = ["8=2I4=2D", "100=", "100=2X20=4I40=1X"]
    for string, tplist in zip(strings, tuplelists):
        cg = Cigar.from_string(string)
        assert cg.pairs == tplist


@pytest.mark.parametrize("l", [0, 10, 1024,  pytest.param(-5, marks=pytest.mark.xfail(raises=ValueError, strict=True))])
@pytest.mark.parametrize("r", [0, 10, 1024,  pytest.param(-5, marks=pytest.mark.xfail(raises=ValueError, strict=True))])
def test_padding(example_cigar, l, r):
    # this fn only tests lengths, but that should be enough
    lref = example_cigar.get_len(ref=True)
    lqry = example_cigar.get_len(ref=False)

    example_cigar.pad(l, r, clip='S')
    # len should stay constant, as clips are ignored by this
    assert example_cigar.get_len(ref=True) == lref
    assert example_cigar.get_len(ref=False) == lqry
    # but start and end should have clips corresponding to the arguments
    if l > 0:
        assert example_cigar.pairs[0][1] == 'S' and example_cigar.pairs[0][0] == l
    if r > 0:
        assert example_cigar.pairs[-1][1] == 'S' and example_cigar.pairs[-1][0] == r

    # test unpadding
    example_cigar.unpad()
    assert example_cigar.get_len(ref=True) == lref
    assert example_cigar.get_len(ref=False) == lqry


