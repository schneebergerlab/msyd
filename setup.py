#/usr/bin/python3
from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

setup(ext_modules = cythonize([
            Extension('ingest', ['ingest.pyx']),
            Extension('pansyn', ['pansyn.pyx'])
            ]),
        scripts=["main.py"]
        )
