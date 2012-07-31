#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"
__status__ = "Development"

"""Test suite for the multiple_assign_taxonomy.py module."""

from cogent.util.unit_test import TestCase, main

from taxcompare.multiple_assign_taxonomy import assign_taxonomy_multiple_times

class MultipleAssignTaxonomyTests(TestCase):
    """Tests for the multiple_assign_taxonomy.py module."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        pass

    def test_assign_taxonomy_multiple_times(self):
        """Functions correctly using standard valid input data."""
        pass


if __name__ == "__main__":
    main()
