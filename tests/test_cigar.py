#!/bin/python3

import msyd
import msyd.cigar as cigar
from msyd.cigar import Cigar
import pytest

@pytest.fixture(params=["100=", "10=2X4=", "100=20I4X5D10I40=", "10=30D5=2X2=20I10="])
def example_cigar(request):
    return cigar.cigar_from_string(request.param)

# cigars guaranteed to have more than 100 bp on both ref and qry
@pytest.fixture(params=["100=", "10=2X40=1X20=1X25=2X10=", "100=20I3=4X4=5D10I40=", "10=30D5=4X2=20I10="])
def example_long_cigar(request):
    return cigar.cigar_from_string(request.param)

# cigars guaranteed to have fewer than 10 bp on both ref and qry
@pytest.fixture(params=["3=1X5=", "9=", "3=2I4D2=", "3D2=1X2="])
def example_short_cigar(request):
    return cigar.cigar_from_string(request.param)

#def test_from_strings():
#    tuplelists = [ [ [8, '='], [2, 'I'], [4, '='], [2, 'D']], [[100, '=']], [[100, '='], [2, 'X'], [20, '='], [4, 'I'], [40, '='], [1, 'X']] ]
#    strings = ["8=2I4=2D", "100=", "100=2X20=4I40=1X"]
#    for string, tplist in zip(strings, tuplelists):
#        cg = cigar.cigar_from_string(string)
#        assert cg.pairs == tplist


#@pytest.mark.parametrize("l", [0, 10, 1024])#,  pytest.param(-5, marks=pytest.mark.xfail(raises=ValueError, strict=True))])
#@pytest.mark.parametrize("r", [0, 10, 1024])#,  pytest.param(-5, marks=pytest.mark.xfail(raises=ValueError, strict=True))])
#def test_padding(example_cigar, l, r):
#    # this fn only tests lengths, but that should be enough
#    lref = example_cigar.get_len(ref=True)
#    lqry = example_cigar.get_len(ref=False)
#
#    example_cigar.pad(l, r, clip='S')
#    # len should stay constant, as clips are ignored by this
#    assert example_cigar.get_len(ref=True) == lref
#    assert example_cigar.get_len(ref=False) == lqry
#
#    # commented out as underlying vector is inaccessible in python
#    # but start and end should have clips corresponding to the arguments
#    #if l > 0:
#    #    assert example_cigar.tups.front().t == b'S' and example_cigar.tups.front().n == l
#    #if r > 0:
#    #    assert example_cigar.tups.back().t == b'S' and example_cigar.tups.back().n == r
#
#    # test unpadding
#    example_cigar.unpad()
#    assert example_cigar.get_len(ref=True) == lref
#    assert example_cigar.get_len(ref=False) == lqry

def test_unequal(example_cigar, example_short_cigar):
    assert not example_cigar == example_short_cigar
    assert example_cigar != example_short_cigar

def test_equal(example_cigar):
    assert example_cigar == example_cigar
    assert not example_cigar != example_cigar

def test_getlen(example_cigar):
    for cg, lr, lq in [("100=", 100, 100), ("10=3X10=", 23, 23), ("10=10I10=", 20, 30), ("10=10D10=", 30, 20)]:
        assert cigar.cigar_from_string(cg).get_len(ref=True) == lr
        assert cigar.cigar_from_string(cg).get_len(ref=False) == lq

selectively_roll_cigar = lambda x, y: (x, cigar.cigar_from_string(y))
def test_getrem():
    # always skip 10 from end/start
    skip = 10
    cgs = [cigar.cigar_from_string(x) for x in ["100=", "5D10=1X10=5I", "50=2I1X4D3="]]
    startrems = [selectively_roll_cigar(*x) for x in [(10, "90="), (5, "5=1X10=5I"), (10, "40=2I1X4D3=")]]
    endrems = [selectively_roll_cigar(*x) for x in [(10, "90="), (15, "5D10=1X"), (8, "48=")]]

    for cg, startrem, endrem in zip(cgs, startrems, endrems):
        assert startrem == cg.get_removed(skip)
        assert endrem == cg.get_removed(skip, start=False)

@pytest.mark.parametrize("ref", [True, False])
@pytest.mark.parametrize("start", [True, False])
def test_getrem_length_invariant(example_cigar, ref, start):
    """
    Checks that the cigar length is preserved by get_rem.
    That is, if a cigar has len n, after removing x, the new length should be n-x.
    """
    l = example_cigar.get_len(ref=ref)
    for i in [10, 50, 100, 200]:
        if i < l:
            assert(example_cigar.get_removed(i, ref=ref, start=start)[1].get_len(ref=ref) + i == l)


@pytest.mark.parametrize("ref", [True, False])
def test_getrem_inv(example_cigar, ref):
    # tests that reversing a cigar string and then skipping X is the same as getting X 
    cg_rev = example_cigar.reverse()

    invskip, invrem = cg_rev.get_removed(12, ref=ref)
    skip, rem = example_cigar.get_removed(12, ref=ref, start=False)
    assert invskip == skip
    assert invrem.reverse() == rem
    assert invrem == rem.reverse()

