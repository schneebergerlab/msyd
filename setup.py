#/usr/bin/python3
from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

setup(ext_modules = cythonize([
            Extension('ingest', ['ingest.pyx'], include_dirs=[numpy.get_include()]),
            Extension('pansyn', ['pansyn.pyx']),
            Extension('varcmp', ['varcmp.pyx']),
            Extension('util', ['util.pyx']),
            Extension('cigar', ['cigar.pyx'])
            ]),
        scripts=["main.py"]
        )
