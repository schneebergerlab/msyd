#/usr/bin/python3
from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

setup(
        ext_modules = cythonize('pansr/*.pyx'),#[ingest_ext, pyx_ext]),
        scripts=["main.py"],
        include_dirs=[numpy.get_include()]
        )
