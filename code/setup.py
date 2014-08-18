#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2014--, taxcompre development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

__version__ = "0.0.0-dev"

import os
from setuptools import find_packages, setup
from setuptools.extension import Extension

classes = """
    Development Status :: 1 - Planning
    License :: OSI Approved :: BSD License
    Topic :: Software Development :: Libraries
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Bio-Informatics
    Programming Language :: Python
    Programming Language :: Python :: 2.7
    Operating System :: Unix
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
"""
classifiers = [s.strip() for s in classes.split('\n') if s]

description = ('Code supporing a systematic comparison of methods for '
               'assigning a taxonomic origin to marker gene sequences.')

setup(name='taxcompare',
      version=__version__,
      license='BSD',
      description=description,
      long_description=description,
      author="taxcompare development team",
      author_email="gregcaporaso@gmail.com",
      maintainer="taxcompare development team",
      maintainer_email="gregcaporaso@gmail.com",
      url='http://github.com/gregcaporaso/short-read-tax-assignment',
      test_suite='nose.collector',
      packages=find_packages(),
      install_requires=['scikit-bio == 0.2.0', 'ipython[all] >= 2.0.0',
                        'biom-format == 2.0.1'],
      extras_require={'test': ["nose >= 0.10.1"]},
      classifiers=classifiers,
      package_data={}
      )
