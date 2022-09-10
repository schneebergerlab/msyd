#!/bin/python3

import pandas as pd
import pytest

import pansyri.util as util

from pansyri.pansyn import *
from pansyri.classes.coords import Range, Pansyn
from pansyri.classes.cigar import Cigar


### Tests for the `Pansyn` module

@pytest.fixture
def simple_pansyn_nocg():
    ret = simple_pansyn()
    ret.cigars_dict = None
    return ret

@pytest.fixture
def simple_pansyn():
    return get_simple_pansyn()

def get_simple_pansyn():
    return Pansyn(
            Range('test', 1, 'NaN', 1, 100),
            {'test1': Range('test1', 1, 'NaN', 1, 100), 'test2': Range('test2', 1, 'NaN', 1, 100)},
            {'test1': Cigar.from_string('100='), 'test2': Cigar.from_string('100=')})

def intermediate_pansyn_nocg():
    ret = intermediate_pansyn()
    ret.cigars_dict = None
    return ret

@pytest.fixture
def intermediate_pansyn():
    return get_intermediate_pansyn()

def get_intermediate_pansyn():
    return Pansyn(
            Range('test', 1, 'NaN', 101, 200),
            {'test1': Range('test1', 1, 'NaN', 101, 200), 'test2': Range('test2', 1, 'NaN', 101, 200)},
            {'test1': Cigar.from_string('50=1X49='), 'test2': Cigar.from_string('100=')})

@pytest.fixture
def complex_pansyn_nocg():
    ret = complex_pansyn()
    ret.cigars_dict = None
    return ret

@pytest.fixture
def complex_pansyn():
    return get_complex_pansyn()

def get_complex_pansyn():
    return Pansyn(
            Range('test', 1, 'NaN', 1001, 1500),
            {
                'test1': Range('test1', 1, 'NaN', 1101, 1600),
                'test2': Range('test2', 2, 'NaN', 901, 1400),
                'test3': Range('test3', 1, 'NaN', 1001, 1550),
                'test4': Range('test4', 3, 'NaN', 951, 1430),
                'test5': Range('test5', 3, 'NaN', 951, 1430)
            },
            {
                'test1': Cigar.from_string('500='),
                'test2': Cigar.from_string('90=10I200=1X99=10D100='),
                'test3': Cigar.from_string('100=20I100=20I200=10I100='),
                'test4': Cigar.from_string('80=10D300=10D100='),
                'test5': Cigar.from_string('100=10I19=1X9=20D50=1X19=1X9=1X30=10I20=20D100=1X49=1X49=')                
            })

@pytest.fixture(params=[get_simple_pansyn, get_intermediate_pansyn, get_complex_pansyn])
def all_pansyns(request):
    return request.param()


# rudimentary test for the __add__ function
def test_add(simple_pansyn, complex_pansyn):
    ret = simple_pansyn + complex_pansyn
    assert ret.ref == simple_pansyn.ref
    assert set(ret.ranges_dict.keys()) == set(simple_pansyn.ranges_dict.keys()).union(complex_pansyn.ranges_dict.keys())
    assert set(ret.cigars_dict.keys()) == set(simple_pansyn.cigars_dict.keys()).union(complex_pansyn.cigars_dict.keys())

@pytest.mark.parametrize("start", [0, 10, 40,  pytest.param(-5, marks=pytest.mark.xfail(raises=ValueError, strict=True))])
@pytest.mark.parametrize("end", [0, 10, 40,  pytest.param(-5, marks=pytest.mark.xfail(raises=ValueError, strict=True))])
def test_drop(all_pansyns, start, end):
    drp = all_pansyns.drop(start, end)
    l = len(drp.ref)
    assert l == len(all_pansyns.ref) - start - end

    for org in all_pansyns.get_organisms():
        assert drp.cigars_dict[org].get_len(ref=False) == len(drp.ranges_dict[org])
        assert drp.cigars_dict[org].get_len(ref=True) == l

@pytest.fixture
def overlapping_pansyns():
    return (Pansyn(Range('test', 1, 'NaN', 101, 200),
        {
            'test1': Range('test1', 1, 'NaN', 151, 250),
            'test2': Range('test2', 1, 'NaN', 101, 200)
        },{
            'test1': Cigar.from_string('10I30=1X29=1X9=10D10=1X9='),
            'test2': Cigar.from_string('30=2X28=1X29=1X4=1X4=')
            }
        ), Pansyn(Range('test', 1, 'NaN', 151, 220),
            {
                'test3': Range('test3', 1, 'NaN', 151, 220),
                'test4': Range('test4', 1, 'NaN', 171, 220)
            }, {
                'test3': Cigar.from_string('30=1X29=1X9='),
                'test4': Cigar.from_string('20=2X28=20D')
                }
            ))


@pytest.mark.parametrize("only_core", [False, True])
@pytest.mark.parametrize("cg", [False, True])
def test_find_overlaps_basic(overlapping_pansyns, cg, only_core):
    """
        Tests some basic properties of calc_overlap output for all possible configurations.
    """
    l, r = overlapping_pansyns
    if not cg:
        l.cigars_dict = None
        r.cigars_dict = None
    res = find_overlaps(pd.DataFrame([l]), pd.DataFrame([r]), only_core=only_core)
    res = [row[1] for row in res.itertuples()]

    # check for proper max sequence count
    assert len(res) <= 3
    # check for proper lengths
    maxlen = max(len(l.ref), len(r.ref))
    for pan in res:
        assert maxlen >= len(pan.ref)# >= MIN_SYN_THRESH

    # check that output is sorted
    assert res == sorted(res)
    
    # check that no organisms are accidentally added
    maxorganisms = set(l.get_organisms()).union(set(r.get_organisms()))
    for pan in res:
        assert set(pan.get_organisms()).issubset(maxorganisms)

    if cg: # if has cgs, check that their length is correct
        for pan in res:
            for org in pan.get_organisms():
                assert len(pan.ref) == pan.cigars_dict[org].get_len(ref=True)
                assert len(pan.ranges_dict[org]) == pan.cigars_dict[org].get_len(ref=False)


ov_ov_nocg = Pansyn(Range('test', 1, 'NaN', 151, 200),
        {
            'test1': Range('test1', 1, 'NaN', 201, 250),
            'test2': Range('test2', 1, 'NaN', 151, 200),
            'test3': Range('test3', 1, 'NaN', 151, 200),
            'test4': Range('test4', 1, 'NaN', 171, 200)
            }, None)

ov_ov = Pansyn(Range('test', 1, 'NaN', 151, 200),
        {
            'test1': Range('test1', 1, 'NaN', 211, 250), # need to drop 50 at start
            'test2': Range('test2', 1, 'NaN', 151, 200),
            'test3': Range('test3', 1, 'NaN', 151, 200), # need to drop 20 at end
            'test4': Range('test4', 1, 'NaN', 171, 220)
        },{
            'test1': Cigar.from_string('10=1X9=10D10=1X9='),
            'test2': Cigar.from_string('10=1X29=1X4=1X4='),
            'test3': Cigar.from_string('30=1X19='),
            'test4': Cigar.from_string('20=2X28=')
        })

ov_noov_left_nocg = Pansyn(Range('test', 1, 'NaN', 101, 150),
        {
            'test1': Range('test1', 1, 'NaN', 151, 200),
            'test2': Range('test2', 1, 'NaN', 101, 150)
            }, None)

ov_noov_right_nocg = Pansyn(Range('test', 1, 'NaN', 201, 220),
        {
            'test3': Range('test3', 1, 'NaN', 201, 220)
            # test4 should be filtered out b/c it would be too small
            }, None)

ov_ov_left = Pansyn(Range('test', 1, 'NaN', 101, 150),
        {
            'test1': Range('test1', 1, 'NaN', 151, 210),
            'test2': Range('test2', 1, 'NaN', 101, 150)
        }, {
            'test1': Cigar.from_string('10I30=1X19='),
            'test2': Cigar.from_string('30=2X18=')
            })

ov_ov_right = Pansyn(Range('test', 1, 'NaN', 201, 220),
        {
            'test3': Range('test3', 1, 'NaN', 201, 220)
            # test4 should be dropped
        }, {
            'test3': Cigar.from_string('10=1X9=')
            })

def test_find_overlaps_nocg(overlapping_pansyns):
    """
        Tests the concrete case of the Pansyns generated by overlapping_pansyns without CIGAR strings.
    """
    l, r = overlapping_pansyns
    l.cigars_dict = None
    r.cigars_dict = None
    l = pd.DataFrame([l])
    r = pd.DataFrame([r])
    res = find_overlaps(l, r, only_core=True)
    ov = res.loc[0][0]
    assert ov == ov_ov_nocg

def test_find_overlaps(overlapping_pansyns):
    """
        Tests the concrete case of the Pansyns generated by overlapping_pansyns with CIGAR strings.
    """
    l, r = overlapping_pansyns
    l = pd.DataFrame([l])
    r = pd.DataFrame([r])
    res = find_overlaps(l, r, only_core=True)
    ov = res.loc[0][0]
    assert ov == ov_ov

def test_find_overlaps_cross_nocg(overlapping_pansyns):
    """
        Tests crosssynteny calling in the concrete case of the Pansyns generated by overlapping_pansyns without CIGAR strings.
    """
    l, r = overlapping_pansyns
    l.cigars_dict = None
    r.cigars_dict = None
    l = pd.DataFrame([l])
    r = pd.DataFrame([r])
    res = find_overlaps(l, r, only_core=False)
    res = [row[1] for row in res.itertuples()]

    assert len(res) == 3
    assert res[1] == ov_ov_nocg # middle should be overlap

    assert res[0] == ov_noov_left_nocg
    assert res[2] == ov_noov_right_nocg

def test_find_overlaps_cross(overlapping_pansyns):
    """
    Tests crossynteny calling in the concrete case of the Pansyns generated by overlapping_pansyns.
    """
    l, r = overlapping_pansyns
    l = pd.DataFrame([l])
    r = pd.DataFrame([r])
    res = find_overlaps(l, r, only_core=False)
    res = [row[1] for row in res.itertuples()]

    assert len(res) == 3
    assert res[1] == ov_ov

    assert res[0] == ov_ov_left
    assert res[2] == ov_ov_right
