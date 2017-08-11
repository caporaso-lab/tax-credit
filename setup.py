#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2014--, tax-credit development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import find_packages, setup

setup(
    name='tax-credit',
    version='0.0.0-dev',
    license='BSD-3-Clause',
    packages=find_packages('tax_credit'),
    install_requires=['biom-format', 'pandas', 'statsmodels', 'bokeh',
                      'scipy', 'jupyter', 'scikit-bio', 'seaborn',
		      'scikit-learn', 'joblib'],
    author="Nicholas Bokulich",
    author_email="nbokulich@gmail.com",
    description="Systematic benchmarking of taxonomic classification methods",
    url="https://github.com/caporaso-lab/tax-credit"
)
