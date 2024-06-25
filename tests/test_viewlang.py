#!/bin/python3

from msyd.util import *
import msyd.util as util
from msyd.coords import Range, Pansyn
import pytest


def test_degree_filters(): 
    ps = Pansyn(Range('ref', 'Chr3', 'x', 100, 200), {
        'org1': Range('org1', 'Chr3','x', 110, 210),
        'org2': Range('org2', 'Chr3','x', 150, 250),
        'org3': Range('org3', 'Chr3','x', 90, 190),
        'org4': Range('org4', 'Chr3','x', 100, 200),
        }, None)
    assert(compile_filter_py('deg >= 2')(ps) == True)
    assert(compile_filter_py('degree>=4')(ps) == True)
    assert(compile_filter_py('deg<=4')(ps) == True)
    assert(compile_filter_py('deg<=7')(ps) == True)
    assert(compile_filter_py('deg<=2')(ps) == False)
    assert(compile_filter_py('deg>=7')(ps) == False)

def test_range_filters():
    ps = Pansyn(Range('ref', 'Chr3', 'x', 100, 200), {
        'org1': Range('org1', 'Chr3','x', 110, 210),
        'org2': Range('org2', 'Chr3','x', 150, 250),
        'org3': Range('org3', 'Chr3','x', 90, 190),
        'org4': Range('org4', 'Chr3','x', 100, 200),
        }, None)
    assert(compile_filter_py('in Chr3:x:50:250')(ps) == True)
    assert(compile_filter_py('in Chr3:x:100:200')(ps) == True) # ends should both be inclusive?
    assert(compile_filter_py('in Chr3:x:99:199')(ps) == False)
    assert(compile_filter_py('in Chr3:x:300:1000')(ps) == False)

    assert(compile_filter_py('on Chr3')(ps) == True)
    assert(compile_filter_py('on Chr4')(ps) == False)

def test_orgs_filters():
    ps = Pansyn(Range('ref', 'Chr3', 'x', 100, 200), {
        'org1': Range('org1', 'Chr3','x', 110, 210),
        'org2': Range('org2', 'Chr3','x', 150, 250),
        'org3': Range('org3', 'Chr3','x', 90, 190),
        'org4': Range('org4', 'Chr3','x', 100, 200),
        }, None)
    assert(compile_filter_py('cont org1')(ps) == True)
    assert(compile_filter_py('contains org9')(ps) == False)
    assert(compile_filter_py('cont ref')(ps) == False)

def test_multi_orgs():
    ps = Pansyn(Range('ref', 'Chr3', 'x', 100, 200), {
        'org1': Range('org1', 'Chr3','x', 110, 210),
        'org2': Range('org2', 'Chr3','x', 150, 250),
        'org3': Range('org3', 'Chr3','x', 90, 190),
        'org4': Range('org4', 'Chr3','x', 100, 200),
        }, None)
    assert(compile_filter_py('cont any org1,org99')(ps) == True)
    assert(compile_filter_py('contains any org98, org1')(ps) == True)
    assert(compile_filter_py('contains any org98, org8')(ps) == False)
    assert(compile_filter_py('cont all org1,org2')(ps) == True)
    assert(compile_filter_py('contains all org1, org2')(ps) == True)
    assert(compile_filter_py('cont all org1,org99')(ps) == False)

def test_not():
    ps = Pansyn(Range('ref', 'Chr3', 'x', 100, 200), {
        'org1': Range('org1', 'Chr3','x', 110, 210),
        'org2': Range('org2', 'Chr3','x', 150, 250),
        'org3': Range('org3', 'Chr3','x', 90, 190),
        'org4': Range('org4', 'Chr3','x', 100, 200),
        }, None)
    assert(compile_filter_py('not cont any org1,org99')(ps) == False)
    assert(compile_filter_py('not on Chr4')(ps) == True)

def test_and():
    ps = Pansyn(Range('ref', 'Chr3', 'x', 100, 200), {
        'org1': Range('org1', 'Chr3','x', 110, 210),
        'org2': Range('org2', 'Chr3','x', 150, 250),
        'org3': Range('org3', 'Chr3','x', 90, 190),
        'org4': Range('org4', 'Chr3','x', 100, 200),
        }, None)
    assert(compile_filter_py('(on Chr5)and(contains org9)')(ps) == False)
    assert(compile_filter_py('(on Chr3) &(contains org9)')(ps) == False)
    assert(compile_filter_py('(not on Chr4) and (on Chr3)')(ps) == True)

def test_or():
    ps = Pansyn(Range('ref', 'Chr3', 'x', 100, 200), {
        'org1': Range('org1', 'Chr3','x', 110, 210),
        'org2': Range('org2', 'Chr3','x', 150, 250),
        'org3': Range('org3', 'Chr3','x', 90, 190),
        'org4': Range('org4', 'Chr3','x', 100, 200),
        }, None)
    assert(compile_filter_py('(on Chr5)or(contains org9)')(ps) == False)
    assert(compile_filter_py('(on Chr3) |(contains org9)')(ps) == True)
    assert(compile_filter_py('(not on Chr4) or (on Chr3)')(ps) == True)

def test_xor():
    ps = Pansyn(Range('ref', 'Chr3', 'x', 100, 200), {
        'org1': Range('org1', 'Chr3','x', 110, 210),
        'org2': Range('org2', 'Chr3','x', 150, 250),
        'org3': Range('org3', 'Chr3','x', 90, 190),
        'org4': Range('org4', 'Chr3','x', 100, 200),
        }, None)
    assert(compile_filter_py('(on Chr5)xor(contains org9)')(ps) == False)
    assert(compile_filter_py('(on Chr3) ^(contains org9)')(ps) == True)
    assert(compile_filter_py('(not on Chr4) and (on Chr3)')(ps) == False)
