#/usr/bin/env python3
import os
os.chdir('./src')
# load version number, follows https://stackoverflow.com/a/16084844
exec(open('msyd/version.py').read())

from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy
from Cython.Compiler import Options


setup(name="msyd",
      version=__version__,
      description='Synteny identification across multiple genome assemblies',
      url='https://github.com/schneebergerlab/msyd/',
      license='MIT License',
      license_files=('LICENSE',),
      ext_modules=cythonize([Extension('msyd.classes.cigar', ['msyd/pyxfiles/classes/cigar.pyx']),
                            Extension('msyd.classes.coords', ['msyd/pyxfiles/classes/coords.pyx']),
                            Extension('msyd.classes.vars', ['msyd/pyxfiles/classes/vars.pyx']),
                            Extension('msyd.io', ['msyd/pyxfiles/io.pyx']),
                            Extension('msyd.pansyn', ['msyd/pyxfiles/pansyn.pyx']),
                            Extension('msyd.realignment', ['msyd/pyxfiles/realignment.pyx']),
                            Extension('msyd.util', ['msyd/pyxfiles/util.pyx']),
                            Extension('msyd.imputation', ['msyd/pyxfiles/imputation.pyx']),
                            Extension('msyd.varcmp', ['msyd/pyxfiles/varcmp.pyx'])]),
      packages=["msyd", "msyd.scripts"],
      include_dirs=[numpy.get_include()],
      scripts=['msyd/scripts/msyd'],
      long_description=open('../README.md').read(),
      )
