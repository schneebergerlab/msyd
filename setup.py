#/usr/bin/env python3
import os
os.chdir('./src')

from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy
#from pasy import __version__
from Cython.Compiler import Options


setup(name="pasy",
      version='0.8.0',
      description='Synteny identification across multiple genome assemblies',
      url='https://github.com/schneebergerlab/pasy/',
      license='MIT License',
      license_files=('LICENSE',),
      ext_modules=cythonize([Extension('pasy.classes.cigar', ['pasy/pyxfiles/classes/cigar.pyx']),
                            Extension('pasy.classes.coords', ['pasy/pyxfiles/classes/coords.pyx']),
                            Extension('pasy.classes.vars', ['pasy/pyxfiles/classes/vars.pyx']),
                             Extension('pasy.io', ['pasy/pyxfiles/io.pyx']),
                             Extension('pasy.pansyn', ['pasy/pyxfiles/pansyn.pyx']),
                             Extension('pasy.util', ['pasy/pyxfiles/util.pyx']),
                             Extension('pasy.ordering', ['pasy/pyxfiles/ordering.pyx']),
                             Extension('pasy.imputation', ['pasy/pyxfiles/imputation.pyx']),
                             Extension('pasy.varcmp', ['pasy/pyxfiles/varcmp.pyx'])]),
      packages=["pasy", "pasy.scripts"],
      include_dirs=[numpy.get_include()],
      scripts=['pasy/scripts/pasy'],
      long_description=open('../README.md').read(),
      )
