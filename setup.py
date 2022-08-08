#/usr/bin/env python3
from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy
from pansr import __version__

setup(name="pansr",
      version='{}'.format(__version__),
      description='Synteny and rearrangement identifier between multiple whole-genome assemblies',
      url='https://github.com/schneebergerlab/pansr/',
      license='MIT License',
      license_files=('LICENSE',),
      ext_modules=cythonize([Extension('pansr.cigar', ['pansr/pyxfiles/cigar.pyx']),
                             Extension('pansr.coords', ['pansr/pyxfiles/coords.pyx']),
                             Extension('pansr.ingest', ['pansr/pyxfiles/ingest.pyx']),
                             Extension('pansr.pansyn', ['pansr/pyxfiles/pansyn.pyx']),
                             Extension('pansr.util', ['pansr/pyxfiles/util.pyx']),
                             Extension('pansr.varcmp', ['pansr/pyxfiles/varcmp.pyx'])]),
      packages=["pansr", "pansr.scripts"],
      include_dirs=[numpy.get_include()],
      scripts=['bin/pansyri'],
      # scripts=['bin/pansyri', 'bin/chroder'],
      long_description=open('README.md').read(),
      )
#
#
# setup(
#         ext_modules = cythonize('pansr/*.pyx'),#[ingest_ext, pyx_ext]),
#         scripts=["main.py"],
#         include_dirs=[numpy.get_include()]
#         )
