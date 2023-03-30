#/usr/bin/env python3
import os
os.chdir('./src')

from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy
#from pansyri import __version__
from Cython.Compiler import Options


setup(name="pansyn",
      version='0.8.0',
      description='Synteny identification across multiple genome assemblies',
      url='https://github.com/schneebergerlab/pansyri/',
      license='MIT License',
      license_files=('LICENSE',),
      ext_modules=cythonize([Extension('pansyn.classes.cigar', ['pansyn/pyxfiles/classes/cigar.pyx']),
                            Extension('pansyn.classes.coords', ['pansyn/pyxfiles/classes/coords.pyx']),
                            Extension('pansyn.classes.vars', ['pansyn/pyxfiles/classes/vars.pyx']),
                             Extension('pansyn.io', ['pansyn/pyxfiles/io.pyx']),
                             Extension('pansyn.pansyn', ['pansyn/pyxfiles/pansyn.pyx']),
                             Extension('pansyn.util', ['pansyn/pyxfiles/util.pyx']),
                             Extension('pansyn.ordering', ['pansyn/pyxfiles/ordering.pyx']),
                             Extension('pansyn.imputation', ['pansyn/pyxfiles/imputation.pyx']),
                             Extension('pansyn.varcmp', ['pansyn/pyxfiles/varcmp.pyx'])]),
      packages=["pansyn", "pansyn.scripts"],
      include_dirs=[numpy.get_include()],
      scripts=['pansyn/scripts/pansyn'],
      long_description=open('../README.md').read(),
      )
