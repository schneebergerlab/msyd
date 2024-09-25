#!/bin/python3

import pandas as pd
import pytest
import copy

import msyd.util as util

from msyd.multisyn import *
from msyd.coords import Range
from msyd.cigar import Cigar
import msyd.cigar as cigar


### Tests for the `Multisyn` module

@pytest.fixture
def simple_multisyn_nocg():
    ret = simple_multisyn()
    ret.cigars_dict = None
    return ret

@pytest.fixture
def simple_multisyn():
    return get_simple_multisyn()

def get_simple_multisyn():
    return Multisyn(
            Range('test', 'chr1', 'NaN', 'chr1', 100),
            {'test1': Range('test1', 'chr1', 'NaN', 'chr1', 100), 'test2': Range('test2', 'chr1', 'NaN', 'chr1', 100)},
            {'test1': cigar.cigar_from_string('100='), 'test2': cigar.cigar_from_string('100=')})

def intermediate_multisyn_nocg():
    ret = intermediate_multisyn()
    ret.cigars_dict = None
    return ret

@pytest.fixture
def intermediate_multisyn():
    return get_intermediate_multisyn()

def get_intermediate_multisyn():
    return Multisyn(
            Range('test', 'chr1', 'NaN', 101, 200),
            {'test1': Range('test1', 'chr1', 'NaN', 101, 200), 'test2': Range('test2', 'chr1', 'NaN', 101, 200)},
            {'test1': cigar.cigar_from_string('50=1X49='), 'test2': cigar.cigar_from_string('100=')})

@pytest.fixture
def complex_multisyn_nocg():
    ret = complex_multisyn()
    ret.cigars_dict = None
    return ret

@pytest.fixture
def complex_multisyn():
    return get_complex_multisyn()

def get_complex_multisyn():
    return Multisyn(
            Range('test', 'chr1', 'NaN', 1001, 1500),
            {
                'test1': Range('test1', 'chr1', 'NaN', 1101, 1600),
                'test2': Range('test2', 'chr2', 'NaN', 901, 1400),
                'test3': Range('test3', 'chr1', 'NaN', 1001, 1550),
                'test4': Range('test4', 'chr3', 'NaN', 951, 1430),
                'test5': Range('test5', 'chr3', 'NaN', 951, 1430)
            },
            {
                'test1': cigar.cigar_from_string('500='),
                'test2': cigar.cigar_from_string('90=10I200=1X99=10D100='),
                'test3': cigar.cigar_from_string('100=20I100=20I200=10I100='),
                'test4': cigar.cigar_from_string('80=10D300=10D100='),
                'test5': cigar.cigar_from_string('100=10I19=1X9=20D50=1X19=1X9=1X30=10I20=20D100=1X49=1X49=')                
            })

@pytest.fixture(params=[get_simple_multisyn, get_intermediate_multisyn, get_complex_multisyn])
def all_multisyns(request):
    return request.param()


# rudimentary test for the __add__ function
def test_add(simple_multisyn, complex_multisyn):
    ret = simple_multisyn + complex_multisyn
    assert ret.ref == simple_multisyn.ref
    assert set(ret.ranges_dict.keys()) == set(simple_multisyn.ranges_dict.keys()).union(complex_multisyn.ranges_dict.keys())
    assert set(ret.cigars_dict.keys()) == set(simple_multisyn.cigars_dict.keys()).union(complex_multisyn.cigars_dict.keys())

@pytest.mark.parametrize("start", [0, 10, 40,  pytest.param(-5, marks=pytest.mark.xfail(raises=ValueError, strict=True))])
@pytest.mark.parametrize("end", [0, 10, 40,  pytest.param(-5, marks=pytest.mark.xfail(raises=ValueError, strict=True))])
def test_drop(all_multisyns, start, end):
    drp = all_multisyns.drop(start, end)
    l = len(drp.ref)
    assert l == len(all_multisyns.ref) - start - end

    for org in all_multisyns.get_organisms():
        assert drp.cigars_dict[org].get_len(ref=False) == len(drp.ranges_dict[org])
        assert drp.cigars_dict[org].get_len(ref=True) == l

@pytest.mark.parametrize("start", [0, 10, 40,  pytest.param(-5, marks=pytest.mark.xfail(raises=ValueError, strict=True))])
@pytest.mark.parametrize("end", [0, 10, 40,  pytest.param(-5, marks=pytest.mark.xfail(raises=ValueError, strict=True))])
def test_dropinplace(all_multisyns, start, end):

    multisyn = copy.copy(all_multisyns)
    #drpctl = multisyn.drop(start, end)
    multisyn.drop_inplace(start, end)


    l = len(multisyn.ref)
    assert l == len(all_multisyns.ref) - start - end

    for org in all_multisyns.get_organisms():
        assert multisyn.cigars_dict[org].get_len(ref=False) == len(multisyns.ranges_dict[org])
        assert multisyn.cigars_dict[org].get_len(ref=True) == l

@pytest.mark.parametrize("start", [0, 10, 40,  pytest.param(-5, marks=pytest.mark.xfail(raises=ValueError, strict=True))])
@pytest.mark.parametrize("end", [0, 10, 40,  pytest.param(-5, marks=pytest.mark.xfail(raises=ValueError, strict=True))])
def test_droponorg(all_multisyns, start, end):
    for org in all_multisyns.get_organisms():
        drp = multisyn.drop_on_org(start, end, org)
        l = len(drp.ranges_dict[org])
        assert l == len(all_multisyns.ranges_dict[org]) - start - end

        for org in drp.get_organisms():
            assert drp.cigars_dict[org].get_len(ref=False) == len(drp.ranges_dict[org])
            assert drp.cigars_dict[org].get_len(ref=True) == l

@pytest.fixture
def overlapping_multisyns():
    return (Multisyn(Range('test', 'chr1', 'NaN', 101, 200),
        {
            'test1': Range('test1', 'chr1', 'NaN', 151, 250),
            'test2': Range('test2', 'chr1', 'NaN', 101, 200)
        },{
            'test1': cigar.cigar_from_string('10I30=1X29=1X9=10D10=1X9='),
            'test2': cigar.cigar_from_string('30=2X28=1X29=1X4=1X4=')
            }
        ), Multisyn(Range('test', 'chr1', 'NaN', 151, 220),
            {
                'test3': Range('test3', 'chr1', 'NaN', 151, 220),
                'test4': Range('test4', 'chr1', 'NaN', 171, 220)
            }, {
                'test3': cigar.cigar_from_string('30=1X29=1X9='),
                'test4': cigar.cigar_from_string('20=2X28=20D')
                }
            ))


@pytest.mark.parametrize("only_core", [False, True])
@pytest.mark.parametrize("cg", [False, True])
def test_find_overlaps_basic(overlapping_multisyns, cg, only_core):
    """
        Tests some basic properties of calc_overlap output for all possible configurations.
    """
    l, r = overlapping_multisyns
    if not cg:
        l.cigars_dict = None
        r.cigars_dict = None
    res = find_overlaps(pd.DataFrame([l]), pd.DataFrame([r]), only_core=only_core)
    res = [row[1] for row in res.itertuples()]

    # check for proper max sequence count
    assert len(res) <= 3
    # check for proper lengths
    maxlen = max(len(l.ref), len(r.ref))
    for multi in res:
        assert maxlen >= len(multi.ref)# >= MIN_SYN_THRESH

    # check that output is sorted
    assert res == sorted(res)
    
    # check that no organisms are accidentally added
    maxorganisms = set(l.get_organisms()).union(set(r.get_organisms()))
    for multi in res:
        assert set(multi.get_organisms()).issubset(maxorganisms)

    if cg: # if has cgs, check that their length is correct
        for multi in res:
            for org in multi.get_organisms():
                assert len(multi.ref) == multi.cigars_dict[org].get_len(ref=True)
                assert len(multi.ranges_dict[org]) == multi.cigars_dict[org].get_len(ref=False)


ov_ov_nocg = Multisyn(Range('test', 'chr1', 'NaN', 151, 200),
        {
            'test1': Range('test1', 'chr1', 'NaN', 201, 250),
            'test2': Range('test2', 'chr1', 'NaN', 151, 200),
            'test3': Range('test3', 'chr1', 'NaN', 151, 200),
            'test4': Range('test4', 'chr1', 'NaN', 171, 200)
            }, None)

ov_ov = Multisyn(Range('test','chr1', 'NaN', 151, 200),
        {
            'test1': Range('test1', 'chr1', 'NaN', 211, 250), # need to drop 50 at start
            'test2': Range('test2', 'chr1', 'NaN', 151, 200),
            'test3': Range('test3', 'chr1', 'NaN', 151, 200), # need to drop 20 at end
            'test4': Range('test4', 'chr1', 'NaN', 171, 220)
        },{
            'test1': cigar.cigar_from_string('10=1X9=10D10=1X9='),
            'test2': cigar.cigar_from_string('10=1X29=1X4=1X4='),
            'test3': cigar.cigar_from_string('30=1X19='),
            'test4': cigar.cigar_from_string('20=2X28=')
        })

ov_noov_left_nocg = Multisyn(Range('test', 'chr1', 'NaN', 101, 150),
        {
            'test1': Range('test1', 'chr1', 'NaN', 151, 200),
            'test2': Range('test2', 'chr1', 'NaN', 101, 150)
            }, None)

ov_noov_right_nocg = Multisyn(Range('test', 'chr1', 'NaN', 201, 220),
        {
            'test3': Range('test3', 'chr1', 'NaN', 201, 220)
            # test4 should be filtered out b/c it would be too small
            }, None)

ov_ov_left = Multisyn(Range('test', 'chr1', 'NaN', 101, 150),
        {
            'test1': Range('test1', 'chr1', 'NaN', 151, 210),
            'test2': Range('test2', 'chr1', 'NaN', 101, 150)
        }, {
            'test1': cigar.cigar_from_string('10I30=1X19='),
            'test2': cigar.cigar_from_string('30=2X18=')
            })

ov_ov_right = Multisyn(Range('test', 'chr1', 'NaN', 201, 220),
        {
            'test3': Range('test3', 'chr1', 'NaN', 201, 220)
            # test4 should be dropped
        }, {
            'test3': cigar.cigar_from_string('10=1X9=')
            })

def test_find_overlaps_nocg(overlapping_multisyns):
    """
        Tests the concrete case of the Multisyns generated by overlapping_multisyns without CIGAR strings.
    """
    l, r = overlapping_multisyns
    l.cigars_dict = None
    r.cigars_dict = None
    l = pd.DataFrame([l])
    r = pd.DataFrame([r])
    res = find_overlaps(l, r, only_core=True)
    ov = res.loc[0][0]
    assert ov == ov_ov_nocg

def test_find_overlaps(overlapping_multisyns):
    """
        Tests the concrete case of the Multisyns generated by overlapping_multisyns with CIGAR strings.
    """
    l, r = overlapping_multisyns
    l = pd.DataFrame([l])
    r = pd.DataFrame([r])
    res = find_overlaps(l, r, only_core=True)
    ov = res.loc[0][0]
    assert ov == ov_ov

def test_find_overlaps_cross_nocg(overlapping_multisyns):
    """
        Tests crosssynteny calling in the concrete case of the Multisyns generated by overlapping_multisyns without CIGAR strings.
    """
    l, r = overlapping_multisyns
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

def test_find_overlaps_cross(overlapping_multisyns):
    """
    Tests crossynteny calling in the concrete case of the Multisyns generated by overlapping_multisyns.
    """
    l, r = overlapping_multisyns
    l = pd.DataFrame([l])
    r = pd.DataFrame([r])
    res = find_overlaps(l, r, only_core=False)
    res = [row[1] for row in res.itertuples()]

    assert len(res) == 3
    assert res[1] == ov_ov

    assert res[0] == ov_ov_left
    assert res[2] == ov_ov_right
