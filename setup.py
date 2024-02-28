#!/usr/bin/env python3

import numpy
from setuptools import setup, Extension
from Cython.Build import cythonize
from Cython.Compiler import Options

import glob

setup(name="msyd",
      description='Synteny identification across multiple genome assemblies',
      url='https://github.com/schneebergerlab/msyd/',
      license='MIT License',
      ext_modules=cythonize([
          Extension(f"msyd.{name.split('/')[-1].split('.')[0]}", [name])
          for name in glob.iglob('msyd/pyxfiles/*.pyx')
          ]),
      packages=["msyd", "msyd.scripts"],
      include_dirs=[numpy.get_include()],
      entry_points={"console_scripts": ["msyd=msyd:main"]},
      long_description=open('./README.md').read(),
      )