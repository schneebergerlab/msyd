#/usr/bin/env python3
from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy
from pansyri import __version__

setup(name="pansyri",
      version='{}'.format(__version__),
      description='Synteny and rearrangement identifier between multiple whole-genome assemblies',
      url='https://github.com/schneebergerlab/pansyri/',
      license='MIT License',
      license_files=('LICENSE',),
      ext_modules=cythonize([Extension('pansyri.classes.cigar', ['pansyri/pyxfiles/classes/cigar.pyx']),
                            Extension('pansyri.classes.coords', ['pansyri/pyxfiles/classes/coords.pyx']),
                             Extension('pansyri.ingest', ['pansyri/pyxfiles/ingest.pyx']),
                             Extension('pansyri.pansyn', ['pansyri/pyxfiles/pansyn.pyx']),
                             Extension('pansyri.util', ['pansyri/pyxfiles/util.pyx']),
                             Extension('pansyri.ordering', ['pansyri/pyxfiles/ordering.pyx']),
                             Extension('pansyri.imputation', ['pansyri/pyxfiles/imputation.pyx']),
                             Extension('pansyri.varcmp', ['pansyri/pyxfiles/varcmp.pyx'])]),
      packages=["pansyri", "pansyri.scripts"],
      include_dirs=[numpy.get_include()],
      scripts=['bin/pansyri'],
      # scripts=['bin/pansyri', 'bin/chroder'],
      long_description=open('README.md').read(),
      )
#
#
# setup(
#         ext_modules = cythonize('pansyri/*.pyx'),#[ingest_ext, pyx_ext]),
#         scripts=["main.py"],
#         include_dirs=[numpy.get_include()]
#         )
