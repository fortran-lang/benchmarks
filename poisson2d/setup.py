#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 16:09:12 2017

@author: arc
"""

from distutils.core import setup
from Cython.Build import cythonize
import numpy
setup(
  name = 'yoo',
  ext_modules = cythonize("poisson.pyx", annotate=True, language_level=3),
  include_dirs=[numpy.get_include()],
)
