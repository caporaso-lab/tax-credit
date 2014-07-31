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

from unittest import TestCase, main
import json
from StringIO import StringIO
import numpy as np

from biom import Table
from taxcompare.eval_framework import (compute_prf, filter_table,
                                       get_observed_observation_ids,
                                       get_actual_and_expected_vectors)

class EvalFrameworkTests(TestCase):

    def setUp(self):
        self.table1 = Table.from_json(json.load(StringIO(table1)))
        self.table2 = Table.from_json(json.load(StringIO(table2)))
        self.table3 = Table.from_json(json.load(StringIO(table3)))

    def test_filter_table(self):
        # prior to filtering there are observations with count less than 10
        self.assertTrue(np.array([e.sum() < 10
            for e in self.table3.iter_data(axis='observation')]).any())
        filtered_table = filter_table(self.table3, min_count=10,
            taxonomy_level=0)
        # after filtering there are no observations with count less than 10
        self.assertFalse(np.array([e.sum() < 10
            for e in filtered_table.iter_data(axis='observation')]).any())

        self.assertTrue(np.array([e.sum() < 100
            for e in self.table3.iter_data(axis='observation')]).any())
        filtered_table = filter_table(self.table3, min_count=100,
            taxonomy_level=0)
        self.assertFalse(np.array([e.sum() < 100
            for e in filtered_table.iter_data(axis='observation')]).any())

        # prior to filtering, there are taxonomies with fewer than 4 levels
        md_levels = [len(md['taxonomy']) < 4
            for _, _, md in self.table3.iter(axis='observation')]
        self.assertTrue(np.array(md_levels).any())
        filtered_table = filter_table(self.table3, min_count=0,
            taxonomy_level=4)
        # after filtering, there are no taxonomies with fewer than 4 levels
        md_levels = [len(md['taxonomy']) < 4
            for _, _, md in filtered_table.iter(axis='observation')]
        self.assertFalse(np.array(md_levels).any())

        md_levels = [len(md['taxonomy']) < 5
            for _, _, md in self.table3.iter(axis='observation')]
        self.assertTrue(np.array(md_levels).any())
        filtered_table = filter_table(self.table3, min_count=0,
            taxonomy_level=5)
        md_levels = [len(md['taxonomy']) < 5
            for _, _, md in filtered_table.iter(axis='observation')]
        self.assertFalse(np.array(md_levels).any())

        md_levels = [len(md['taxonomy']) < 6
            for _, _, md in self.table3.iter(axis='observation')]
        self.assertTrue(np.array(md_levels).any())
        filtered_table = filter_table(self.table3, min_count=0,
            taxonomy_level=6)
        md_levels = [len(md['taxonomy']) < 6
            for _, _, md in filtered_table.iter(axis='observation')]
        self.assertFalse(np.array(md_levels).any())

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

    def test_get_actual_and_expected_vectors(self):
        actual_actual_vector, actual_expected_vector = \
         get_actual_and_expected_vectors(self.table1, self.table2)

        # We don't care about the order of the vector, just the
        # pairing of the values across the two vectors.
        actual = zip(actual_actual_vector, actual_expected_vector)
        actual.sort()

        expected_actual_vector = [1.0, 2.0, 9.0]
        expected_expected_vector = [1.0, 0.0, 9.0]
        expected = zip(expected_actual_vector,expected_expected_vector)
        expected.sort()

        self.assertEqual(actual,expected)

        ## look at samples other than the first in each table
        actual_actual_vector, actual_expected_vector = \
         get_actual_and_expected_vectors(self.table1,
                                         self.table2,
                                         actual_sample_id='s2',
                                         expected_sample_id='s4')

        # We don't care about the order of the vector, just the
        # pairing of the values across the two vectors.
        actual = zip(actual_actual_vector, actual_expected_vector)
        actual.sort()

        # o2 is not observed, so its value shouldn't show up in the vectors
        expected_actual_vector = [2.0, 10.0]
        expected_expected_vector = [0.001, 0.0]
        expected = zip(expected_actual_vector,expected_expected_vector)
        expected.sort()

        self.assertEqual(actual,expected)

table1 = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "greg","date": "2013-08-22T13:10:23.907145","matrix_type": "sparse","matrix_element_type": "float","shape": [3, 4],"data": [[0,0,1.0],[0,1,2.0],[0,2,3.0],[0,3,4.0],[1,0,2.0],[1,1,0.0],[1,2,7.0],[1,3,8.0],[2,0,9.0],[2,1,10.0],[2,2,11.0],[2,3,12.0]],"rows": [{"id": "o1", "metadata": {"domain": "Archaea"}},{"id": "o2", "metadata": {"domain": "Bacteria"}},{"id": "o3", "metadata": {"domain": "Bacteria"}}],"columns": [{"id": "s1", "metadata": {"country": "Peru", "pH": 4.2}},{"id": "s2", "metadata": {"country": "Peru", "pH": 5.2}},{"id": "s3", "metadata": {"country": "Peru", "pH": 5.0}},{"id": "s4", "metadata": {"country": "Peru", "pH": 4.9}}]}"""

## table 1
#OTU ID	s1	s2	s3	s4
# o1    1.0 2.0 3.0 4.0
# o2    2.0 0.0 7.0 8.0
# o3    9.0 10.0    11.0    12.0


table2 = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "greg","date": "2013-08-22T13:19:35.281188","matrix_type": "sparse","matrix_element_type": "float","shape": [2, 4],"data": [[0,0,1.0],[0,1,2.0],[0,2,3.0],[0,3,0.001],[1,0,9.0],[1,1,10.0],[1,2,11.0],[1,3,0.0]],"rows": [{"id": "o1", "metadata": {"domain": "Archaea"}},{"id": "o3", "metadata": {"domain": "Bacteria"}}],"columns": [{"id": "s1", "metadata": {"country": "Peru", "pH": 4.2}},{"id": "s2", "metadata": {"country": "Peru", "pH": 5.2}},{"id": "s3", "metadata": {"country": "Peru", "pH": 5.0}},{"id": "s4", "metadata": {"country": "Peru", "pH": 4.9}}]}"""

## table 2
#OTU ID	s1	s2	s3	s4
# o1    1.0 2.0 3.0 0.001
# o3    9.0 10.0    11.0    0.0

table3 = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "BIOM-Format 1.1.2","date": "2013-06-13T09:41:43.709874","matrix_type": "sparse","matrix_element_type": "float","shape": [70, 4],"data": [[0,0,1.0],[0,1,1.0],[1,0,1.0],[1,1,2.0],[1,2,2.0],[1,3,1.0],[2,0,22.0],[2,1,44.0],[2,2,19.0],[2,3,26.0],[3,0,937.0],[3,1,1815.0],[3,2,923.0],[3,3,775.0],[4,0,1.0],[4,1,1.0],[4,2,3.0],[4,3,1.0],[5,0,130.0],[5,1,229.0],[5,2,122.0],[5,3,69.0],[6,2,1.0],[6,3,2.0],[7,0,52.0],[7,1,80.0],[7,2,5.0],[7,3,2.0],[8,1,2.0],[9,0,3.0],[9,1,7.0],[9,2,4.0],[9,3,2.0],[10,1,1.0],[10,3,1.0],[11,0,6.0],[11,1,9.0],[11,2,4.0],[11,3,5.0],[12,1,1.0],[12,2,1.0],[12,3,2.0],[13,0,1.0],[13,2,1.0],[14,1,2.0],[15,0,1.0],[15,3,3.0],[16,3,2.0],[17,1,4.0],[18,0,1.0],[18,3,1.0],[19,0,1.0],[19,1,1.0],[19,3,1.0],[20,0,5.0],[20,1,13.0],[21,0,2.0],[21,1,3.0],[21,2,2.0],[21,3,1.0],[22,0,1.0],[22,1,2.0],[23,0,2.0],[23,1,2.0],[23,2,2.0],[23,3,1.0],[24,0,1.0],[24,1,1.0],[25,1,2.0],[25,3,1.0],[26,0,17.0],[26,1,18.0],[26,2,69.0],[26,3,64.0],[27,1,1.0],[27,3,2.0],[28,0,20.0],[28,1,29.0],[28,2,133.0],[28,3,104.0],[29,0,2.0],[29,1,5.0],[29,2,2.0],[29,3,3.0],[30,0,31.0],[30,1,48.0],[30,2,10.0],[30,3,15.0],[31,0,1.0],[31,1,2.0],[31,2,15.0],[31,3,12.0],[32,0,1.0],[32,1,1.0],[33,0,94.0],[33,1,150.0],[33,2,63.0],[33,3,39.0],[34,1,1.0],[34,2,1.0],[35,1,4.0],[36,0,1.0],[36,1,1.0],[37,1,1.0],[37,3,1.0],[38,0,1.0],[38,1,1.0],[39,0,22.0],[39,1,44.0],[39,2,1.0],[40,0,4.0],[40,1,7.0],[41,0,1.0],[41,1,2.0],[41,2,3.0],[42,0,198.0],[42,1,374.0],[42,2,181.0],[42,3,167.0],[43,0,192.0],[43,1,338.0],[43,2,5.0],[43,3,17.0],[44,0,1.0],[44,1,1.0],[45,0,1.0],[45,1,1.0],[45,3,1.0],[46,0,1.0],[46,1,1.0],[46,3,4.0],[47,0,2.0],[47,1,3.0],[47,2,1.0],[47,3,3.0],[48,1,1.0],[48,2,1.0],[49,0,2.0],[49,1,1.0],[50,0,14.0],[50,1,19.0],[50,2,6.0],[50,3,8.0],[51,0,27.0],[51,1,55.0],[51,2,1.0],[52,1,1.0],[52,2,1.0],[53,2,2.0],[54,0,9.0],[54,1,27.0],[54,2,14.0],[54,3,11.0],[55,1,1.0],[55,3,1.0],[56,0,8.0],[56,1,9.0],[56,2,2.0],[56,3,4.0],[57,0,1.0],[57,1,1.0],[57,2,1.0],[57,3,1.0],[58,0,3.0],[58,1,1.0],[58,2,1.0],[58,3,1.0],[59,1,2.0],[60,0,3.0],[60,2,1.0],[61,0,91.0],[61,1,160.0],[61,2,4.0],[61,3,3.0],[62,0,1.0],[62,1,1.0],[62,2,1.0],[62,3,2.0],[63,0,3.0],[63,1,1.0],[64,0,1.0],[64,1,1.0],[64,2,2.0],[64,3,1.0],[65,2,1.0],[65,3,1.0],[66,1,2.0],[66,2,2.0],[66,3,2.0],[67,2,1.0],[67,3,1.0],[68,0,1.0],[68,1,2.0],[69,0,1.0],[69,1,1.0]],"rows": [{"id": "269901", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Pseudomonadales", "f__Pseudomonadaceae"]}},{"id": "4130483", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "137056", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "1995363", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Bacillales", "f__Staphylococcaceae", "g__Staphylococcus"]}},{"id": "939252", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Bacillales", "f__Staphylococcaceae", "g__Staphylococcus"]}},{"id": "4380971", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Clostridiaceae", "g__Clostridium", "s__butyricum"]}},{"id": "1081058", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Bacillales", "f__Staphylococcaceae", "g__Staphylococcus"]}},{"id": "4440404", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Betaproteobacteria", "o__Neisseriales", "f__Neisseriaceae", "g__Neisseria"]}},{"id": "984924", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "953758", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Bacillales", "f__Staphylococcaceae", "g__Staphylococcus"]}},{"id": "4316928", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "152001", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Bacillales"]}},{"id": "227083", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "4445673", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Clostridiales", "f__Clostridiaceae", "g__Clostridium", "s__perfringens"]}},{"id": "138389", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "4427114", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "153046", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Bacillales", "f__Staphylococcaceae", "g__Staphylococcus"]}},{"id": "1059655", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Lactobacillales", "f__Streptococcaceae", "g__Streptococcus"]}},{"id": "1550056", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "979261", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "12574", "metadata": {"taxonomy": ["k__Bacteria", "p__Actinobacteria", "c__Actinobacteria", "o__Actinomycetales", "f__Actinomycetaceae", "g__Actinomyces", "s__"]}},{"id": "368134", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Bacillales", "f__Staphylococcaceae"]}},{"id": "1039016", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "996487", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Bacillales", "f__Staphylococcaceae", "g__Staphylococcus"]}},{"id": "1069592", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "1112200", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "4297222", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes"]}},{"id": "923151", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "532163", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Alphaproteobacteria", "o__Rhodobacterales", "f__Rhodobacteraceae"]}},{"id": "928538", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Bacillales", "f__Staphylococcaceae", "g__Staphylococcus"]}},{"id": "1891556", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "114510", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Enterobacteriales", "f__Enterobacteriaceae"]}},{"id": "158047", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "242070", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Pseudomonadales", "f__Pseudomonadaceae"]}},{"id": "149265", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "919490", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "164413", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "767863", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "113773", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Lactobacillales"]}},{"id": "128604", "metadata": {"taxonomy": ["k__Bacteria", "p__[Thermi]", "c__Deinococci", "o__Deinococcales", "f__Deinococcaceae", "g__Deinococcus", "s__"]}},{"id": "99882", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Lactobacillales"]}},{"id": "519673", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "630141", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "219151", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Pseudomonadales", "f__Moraxellaceae", "g__Acinetobacter"]}},{"id": "977188", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "1121111", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "894774", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "441155", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Bacillales", "f__Staphylococcaceae", "g__Staphylococcus"]}},{"id": "1059977", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Lactobacillales", "f__Streptococcaceae", "g__Streptococcus"]}},{"id": "552922", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Pseudomonadales", "f__Moraxellaceae", "g__Acinetobacter"]}},{"id": "2874742", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "1756274", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Lactobacillales", "f__Enterococcaceae"]}},{"id": "4315958", "metadata": {"taxonomy": ["k__Bacteria"]}},{"id": "617833", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "2896107", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Bacillales", "f__Staphylococcaceae"]}},{"id": "4365141", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Lactobacillales", "f__Leuconostocaceae", "g__Leuconostoc"]}},{"id": "356733", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "1067519", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "1068955", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "4438739", "metadata": {"taxonomy": ["k__Bacteria", "p__Actinobacteria", "c__Actinobacteria", "o__Actinomycetales", "f__Propionibacteriaceae", "g__Propionibacterium"]}},{"id": "164612", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "4416988", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "1055132", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "187233", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Lactobacillales", "f__Lactobacillaceae"]}},{"id": "New.CleanUp.ReferenceOTU0", "metadata": {"taxonomy": ["k__Bacteria"]}},{"id": "New.CleanUp.ReferenceOTU2", "metadata": {"taxonomy": ["k__Bacteria", "p__Proteobacteria", "c__Alphaproteobacteria"]}},{"id": "New.CleanUp.ReferenceOTU10", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Bacillales", "f__Staphylococcaceae"]}},{"id": "New.CleanUp.ReferenceOTU27", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}},{"id": "New.CleanUp.ReferenceOTU36", "metadata": {"taxonomy": ["k__Bacteria"]}},{"id": "New.CleanUp.ReferenceOTU39", "metadata": {"taxonomy": ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]}}],"columns": [{"id": "HMPMockV1.1.Even1", "metadata": null},{"id": "HMPMockV1.1.Even2", "metadata": null},{"id": "HMPMockV1.2.Staggered2", "metadata": null},{"id": "HMPMockV1.2.Staggered1", "metadata": null}]}"""

if __name__ == "__main__":
    main()
