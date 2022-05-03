#/usr/bin/python3
from setuptools import setup
from Cython.Build import cythonize
import numpy

setup(ext_modules = cythonize("ingest.pyx"),
        scripts=["pansr.pyx"]
        )
