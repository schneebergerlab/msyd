#!/bin/python3

import pandas as pd
import pytest

import pansyri.util as util

from pansyri.pansyn import *
from pansyri.classes.coords import Range, Pansyn
from pansyri.classes.cigar import Cigar


### Tests for the `Pansyn` class

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


def test_calc_overlap():
    #TODO test calc_overlap

    # generate pansyns here, makes more sense than to use standard ones
    left = Pansyn()
    right = Pansyn()


    




