#!/usr/bin/env python

from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from cogent.util.unit_test import TestCase, main
from biom.parse import parse_biom_table
from taxcompare.eval_framework import (compute_prf,
                                       get_observed_observation_ids)

"""Test suite for the eval_framework.py module."""

class EvalFrameworkTests(TestCase):
    
    def setUp(self):
        self.table1 = parse_biom_table(table1.split('\n'))
        self.table2 = parse_biom_table(table2.split('\n'))
    
    def test_compute_prf_default(self):
        """ p, r and f are computed correctly when defaulting to first sample ids"""
        # default of comparing first sample in each table
        actual = compute_prf(self.table1,self.table2)
        expected = (2./3., 1.0, 0.8)
        self.assertEqual(actual,expected)
        # default of comparing first sample in each table
        actual = compute_prf(self.table2,self.table1)
        expected = (1.0, 2./3., 0.8)
        self.assertEqual(actual,expected)

    def test_compute_prf_alt_sample_ids(self):
        """ p, r and f are computed correctly when using alternative sample ids"""
        # alt sample in table 1
        actual = compute_prf(self.table1,self.table2,actual_sample_id='s2')
        expected = (1.0, 1.0, 1.0)
        self.assertEqual(actual,expected)
        
        # alt sample in table 2
        actual = compute_prf(self.table1,self.table2,expected_sample_id='s4')
        expected = (1./3., 1.0, 0.5)
        self.assertEqual(actual,expected)
        
        # alt sample in tables 1 & 2
        actual = compute_prf(self.table1,self.table2,
                             actual_sample_id='s2',expected_sample_id='s4')
        expected = (0.5, 1.0, 2./3.)
        self.assertEqual(actual,expected)
        

table1 = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "greg","date": "2013-08-22T13:10:23.907145","matrix_type": "sparse","matrix_element_type": "float","shape": [3, 4],"data": [[0,0,1.0],[0,1,2.0],[0,2,3.0],[0,3,4.0],[1,0,2.0],[1,1,0.0],[1,2,7.0],[1,3,8.0],[2,0,9.0],[2,1,10.0],[2,2,11.0],[2,3,12.0]],"rows": [{"id": "o1", "metadata": {"domain": "Archaea"}},{"id": "o2", "metadata": {"domain": "Bacteria"}},{"id": "o3", "metadata": {"domain": "Bacteria"}}],"columns": [{"id": "s1", "metadata": {"country": "Peru", "pH": 4.2}},{"id": "s2", "metadata": {"country": "Peru", "pH": 5.2}},{"id": "s3", "metadata": {"country": "Peru", "pH": 5.0}},{"id": "s4", "metadata": {"country": "Peru", "pH": 4.9}}]}"""

table2 = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "greg","date": "2013-08-22T13:19:35.281188","matrix_type": "sparse","matrix_element_type": "float","shape": [2, 4],"data": [[0,0,1.0],[0,1,2.0],[0,2,3.0],[0,3,0.001],[1,0,9.0],[1,1,10.0],[1,2,11.0],[1,3,0.0]],"rows": [{"id": "o1", "metadata": {"domain": "Archaea"}},{"id": "o3", "metadata": {"domain": "Bacteria"}}],"columns": [{"id": "s1", "metadata": {"country": "Peru", "pH": 4.2}},{"id": "s2", "metadata": {"country": "Peru", "pH": 5.2}},{"id": "s3", "metadata": {"country": "Peru", "pH": 5.0}},{"id": "s4", "metadata": {"country": "Peru", "pH": 4.9}}]}"""

if __name__ == "__main__":
    main()