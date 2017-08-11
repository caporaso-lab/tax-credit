#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2014--, tax-credit development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


from unittest import TestCase, main
import json
import numpy as np
import pandas as pd
from io import StringIO
from os.path import join
from tempfile import mkdtemp
from shutil import rmtree

from biom import Table
from sklearn.metrics import precision_recall_fscore_support

from tax_credit.eval_framework import (compute_taxon_accuracy,
                                       filter_table,
                                       get_sample_to_top_params,
                                       parameter_comparisons,
                                       per_sequence_precision)


class EvalFrameworkTests(TestCase):

    def test_get_sample_to_top_params(self):
        actual = get_sample_to_top_params(self.mock_result_table1, "F-measure")
        self.assertEqual(actual['rdp'][('B1', 'm1')], ['0.2'])
        self.assertEqual(actual['rdp'][('F2', 'm2')], ['0.1'])
        self.assertEqual(actual['rdp'][('F2', 'm3')], ['0', '0.1'])
        self.assertEqual(actual['uclust'][('B1', 'm1')], ['0.51:0.8:3'])
        self.assertEqual(actual['uclust'][('F2', 'm2')], ['0.51:0.8:3'])
        self.assertEqual(actual['uclust'][('F2', 'm3')], ['0.51:0.9:3'])
        self.assertEqual(actual.shape, (3, 2))

    def test_parameter_comparisons(self):
        actual = parameter_comparisons(self.mock_result_table1, "rdp")
        self.assertEqual(actual['F-measure']['0.1'], 2)
        self.assertEqual(actual['F-measure']['0.2'], 1)
        self.assertEqual(actual['F-measure']['0'], 1)
        self.assertEqual(actual['F-measure']['0.3'], 0)
        self.assertEqual(actual['Precision']['0.1'], 2)
        self.assertEqual(actual['Recall']['0.1'], 3)
        self.assertEqual(actual.shape, (6, 5))

        actual = parameter_comparisons(self.mock_result_table1, "uclust")
        self.assertEqual(actual['F-measure']['0.51:0.8:3'], 2)
        self.assertEqual(actual['F-measure']['0.51:0.9:3'], 1)
        self.assertEqual(actual.shape, (2, 5))

    def test_filter_table(self):
        # prior to filtering there are observations with count less than 10
        self.assertTrue(np.array(
            [e.sum() < 10 for e in self.table3.iter_data(
             axis='observation')]).any())
        filtered_table = filter_table(
            self.table3, min_count=10, taxonomy_level=0)
        # after filtering there are no observations with count less than 10
        self.assertFalse(np.array(
            [e.sum() < 10 for e in filtered_table.iter_data(
             axis='observation')]).any())
        # but some observations are still present
        self.assertTrue(filtered_table.shape[0] > 0)

        self.assertTrue(np.array(
            [e.sum() < 100 for e in self.table3.iter_data(
             axis='observation')]).any())
        filtered_table = filter_table(
            self.table3, min_count=100, taxonomy_level=0)
        self.assertFalse(np.array(
            [e.sum() < 100 for e in filtered_table.iter_data(
             axis='observation')]).any())
        # but some observations are still present
        self.assertTrue(filtered_table.shape[0] > 0)

        # prior to filtering, there are taxonomies with fewer than 4 levels
        md_levels = [len(md['taxonomy']) < 4
                     for _, _, md in self.table3.iter(axis='observation')]
        self.assertTrue(np.array(md_levels).any())
        filtered_table = filter_table(
            self.table3, min_count=0, taxonomy_level=4)
        # after filtering, there are no taxonomies with fewer than 4 levels
        md_levels = [len(md['taxonomy']) < 4
                     for _, _, md in filtered_table.iter(axis='observation')]
        self.assertFalse(np.array(md_levels).any())
        # but some observations are still present
        self.assertTrue(filtered_table.shape[0] > 0)

        md_levels = [len(md['taxonomy']) < 5
                     for _, _, md in self.table3.iter(axis='observation')]
        self.assertTrue(np.array(md_levels).any())
        filtered_table = filter_table(
            self.table3, min_count=0, taxonomy_level=5)
        md_levels = [len(md['taxonomy']) < 5
                     for _, _, md in filtered_table.iter(axis='observation')]
        self.assertFalse(np.array(md_levels).any())
        # but some observations are still present
        self.assertTrue(filtered_table.shape[0] > 0)

        md_levels = [len(md['taxonomy']) < 6
                     for _, _, md in self.table3.iter(axis='observation')]
        self.assertTrue(np.array(md_levels).any())
        filtered_table = filter_table(
            self.table3, min_count=0, taxonomy_level=6)
        md_levels = [len(md['taxonomy']) < 6
                     for _, _, md in filtered_table.iter(axis='observation')]
        self.assertFalse(np.array(md_levels).any())
        # but some observations are still present
        self.assertTrue(filtered_table.shape[0] > 0)

    def test_filter_table_taxa(self):
        """ taxa-based filtering works as expected """
        taxa_to_keep = ["k__Bacteria", "p__Firmicutes", "c__Bacilli"]
        filtered_table = filter_table(self.table3, taxa_to_keep=taxa_to_keep)
        # expected value determined with grep -c c__Bacilli
        self.assertEqual(filtered_table.shape[0], 53)

        taxa_to_keep = ["k__Bacteria", "p__Firmicutes", "c__Bacilli",
                        "o__Bacillales", "f__Staphylococcaceae",
                        "g__Staphylococcus"]
        filtered_table = filter_table(self.table3, taxa_to_keep=taxa_to_keep)
        # expected value determined with grep -c g__Staphylococcus
        self.assertEqual(filtered_table.shape[0], 8)

        taxa_to_keep = ["k__Bacteria"]
        filtered_table = filter_table(self.table3, taxa_to_keep=taxa_to_keep)
        # all observations are retained
        self.assertEqual(filtered_table.shape[0], self.table3.shape[0])

        taxa_to_keep = ["k__Archaea"]
        filtered_table = filter_table(self.table3, taxa_to_keep=taxa_to_keep)
        # no observations are retained
        self.assertEqual(filtered_table.shape[0], 0)

    def test_compute_taxon_accuracy_default(self):
        """ p, r and f compute correctly when default to first sample ids"""
        # default of comparing first sample in each table
        actual = compute_taxon_accuracy(self.table1, self.table2)
        expected = (2./3., 1.0)
        self.assertAlmostEqual(actual, expected)
        # default of comparing first sample in each table
        actual = compute_taxon_accuracy(self.table2, self.table1)
        expected = (1.0, 2./3.)
        self.assertAlmostEqual(actual, expected)

    def test_compute_taxon_accuracy_alt_sample_ids(self):
        """ p, r and f compute correctly when using alternative sample ids"""
        # alt sample in table 1
        actual = compute_taxon_accuracy(
            self.table1, self.table2, actual_sample_id='s2')
        expected = (1.0, 1.0)
        self.assertEqual(actual, expected)

        # alt sample in table 2
        actual = compute_taxon_accuracy(
            self.table1, self.table2, expected_sample_id='s4')
        expected = (1./3., 1.0)
        self.assertAlmostEqual(actual, expected)

        # alt sample in tables 1 & 2
        actual = compute_taxon_accuracy(
            self.table1, self.table2, actual_sample_id='s2',
            expected_sample_id='s4')
        expected = (0.5, 1.0)
        self.assertAlmostEqual(actual, expected)

    def test_per_sequence_precision(self):
        """ p, r, and f on individual expected sequences."""
        exp_taxa = join(self.tmpdir, 'trueish-taxonomies.tsv')
        obs_taxa = join(self.tmpdir, 'taxonomy.tsv')

        exp_list = []
        obs_list = []
        weights = []
        with open(exp_taxa) as exp_in:
            for line in exp_in:
                sid, taxon = line.split()
                if not self.table1.exists(sid, axis='observation'):
                    continue
                exp_list.append(taxon)
                with open(obs_taxa) as obs_in:
                    for oline in obs_in:
                        osid, otaxon = oline.split()
                        if osid == sid:
                            obs_list.append(otaxon)
                            break
                weights.append(self.table1.get_value_by_ids(sid, 's1'))
        for i in range(7):
            p, r, f = per_sequence_precision(
                exp_taxa, obs_taxa, self.table1, 's1', i)
            elist = [';'.join(e.split(';')[:i+1]) for e in exp_list]
            olist = [';'.join(o.split(';')[:i+1]) for o in obs_list]
            ep, er, ef, _ = precision_recall_fscore_support(
                elist, olist, sample_weight=weights, average='micro')
            self.assertAlmostEqual(p, ep)
            self.assertAlmostEqual(r, er)
            self.assertAlmostEqual(f, ef)

    def test_per_sequence_precision_exclude(self):
        """ p, r, and f on individual expected sequences, excluding specific
        taxa.
        """
        exp_taxa = join(self.tmpdir, 'trueish-taxonomies.tsv')
        obs_taxa = join(self.tmpdir, 'taxonomy.tsv')
        etable = self.table1.filter(['o1'], axis='observation', invert=True,
                                    inplace=False)
        for i in range(7):
            p, r, f = per_sequence_precision(
                exp_taxa, obs_taxa, self.table1, 's1', i, exclude=['other'])
            ep, er, ef = per_sequence_precision(
                exp_taxa, obs_taxa, etable, 's1', i)
            self.assertAlmostEqual(p, ep)
            self.assertAlmostEqual(r, er)
            self.assertAlmostEqual(f, ef)

    @classmethod
    def setUpClass(self):
        _table1 = """{"id": "None",
                      "format": "Biological Observation Matrix 1.0.0",
                      "format_url": "http:\/\/biom-format.org",
                      "type": "OTU table",
                      "generated_by": "greg",
                      "date": "2013-08-22T13:10:23.907145",
                      "matrix_type": "sparse",
                      "matrix_element_type": "float",
                      "shape": [
                        3,
                        4
                      ],
                      "data": [
                        [
                          0,
                          0,
                          1
                        ],
                        [
                          0,
                          1,
                          2
                        ],
                        [
                          0,
                          2,
                          3
                        ],
                        [
                          0,
                          3,
                          4
                        ],
                        [
                          1,
                          0,
                          2
                        ],
                        [
                          1,
                          1,
                          0
                        ],
                        [
                          1,
                          2,
                          7
                        ],
                        [
                          1,
                          3,
                          8
                        ],
                        [
                          2,
                          0,
                          9
                        ],
                        [
                          2,
                          1,
                          10
                        ],
                        [
                          2,
                          2,
                          11
                        ],
                        [
                          2,
                          3,
                          12
                        ]
                      ],
                      "rows": [
                        {
                          "id": "o1",
                          "metadata": {
                            "domain": "Archaea"
                          }
                        },
                        {
                          "id": "o2",
                          "metadata": {
                            "domain": "Bacteria"
                          }
                        },
                        {
                          "id": "o3",
                          "metadata": {
                            "domain": "Bacteria"
                          }
                        }
                      ],
                      "columns": [
                        {
                          "id": "s1",
                          "metadata": {
                            "country": "Peru",
                            "pH": 4.2
                          }
                        },
                        {
                          "id": "s2",
                          "metadata": {
                            "country": "Peru",
                            "pH": 5.2
                          }
                        },
                        {
                          "id": "s3",
                          "metadata": {
                            "country": "Peru",
                            "pH": 5
                          }
                        },
                        {
                          "id": "s4",
                          "metadata": {
                            "country": "Peru",
                            "pH": 4.9
                          }
                        }
                      ]
                    }"""
        # table 1
        # OTU ID	s1	s2	s3	s4
        # o1    1.0 2.0 3.0 4.0
        # o2    2.0 0.0 7.0 8.0
        # o3    9.0 10.0    11.0    12.0

        _table2 = """{"id": "None",
                      "format": "Biological Observation Matrix 1.0.0",
                      "format_url": "http:\/\/biom-format.org",
                      "type": "OTU table",
                      "generated_by": "greg",
                      "date": "2013-08-22T13:19:35.281188",
                      "matrix_type": "sparse",
                      "matrix_element_type": "float",
                      "shape": [
                        2,
                        4
                      ],
                      "data": [
                        [
                          0,
                          0,
                          1
                        ],
                        [
                          0,
                          1,
                          2
                        ],
                        [
                          0,
                          2,
                          3
                        ],
                        [
                          0,
                          3,
                          0.001
                        ],
                        [
                          1,
                          0,
                          9
                        ],
                        [
                          1,
                          1,
                          10
                        ],
                        [
                          1,
                          2,
                          11
                        ],
                        [
                          1,
                          3,
                          0
                        ]
                      ],
                      "rows": [
                        {
                          "id": "o1",
                          "metadata": {
                            "domain": "Archaea"
                          }
                        },
                        {
                          "id": "o3",
                          "metadata": {
                            "domain": "Bacteria"
                          }
                        }
                      ],
                      "columns": [
                        {
                          "id": "s1",
                          "metadata": {
                            "country": "Peru",
                            "pH": 4.2
                          }
                        },
                        {
                          "id": "s2",
                          "metadata": {
                            "country": "Peru",
                            "pH": 5.2
                          }
                        },
                        {
                          "id": "s3",
                          "metadata": {
                            "country": "Peru",
                            "pH": 5
                          }
                        },
                        {
                          "id": "s4",
                          "metadata": {
                            "country": "Peru",
                            "pH": 4.9
                          }
                        }
                      ]
                    }"""

        # table 2
        # OTU ID	s1	s2	s3	s4
        # o1    1.0 2.0 3.0 0.001
        # o3    9.0 10.0    11.0    0.0

        _table3 = """{"id": "None",
                      "format": "Biological Observation Matrix 1.0.0",
                      "format_url": "http:\/\/biom-format.org",
                      "type": "OTU table",
                      "generated_by": "BIOM-Format 1.1.2",
                      "date": "2013-06-13T09:41:43.709874",
                      "matrix_type": "sparse",
                      "matrix_element_type": "float",
                      "shape": [
                        70,
                        4
                      ],
                      "data": [
                        [
                          0,
                          0,
                          1
                        ],
                        [
                          0,
                          1,
                          1
                        ],
                        [
                          1,
                          0,
                          1
                        ],
                        [
                          1,
                          1,
                          2
                        ],
                        [
                          1,
                          2,
                          2
                        ],
                        [
                          1,
                          3,
                          1
                        ],
                        [
                          2,
                          0,
                          22
                        ],
                        [
                          2,
                          1,
                          44
                        ],
                        [
                          2,
                          2,
                          19
                        ],
                        [
                          2,
                          3,
                          26
                        ],
                        [
                          3,
                          0,
                          937
                        ],
                        [
                          3,
                          1,
                          1815
                        ],
                        [
                          3,
                          2,
                          923
                        ],
                        [
                          3,
                          3,
                          775
                        ],
                        [
                          4,
                          0,
                          1
                        ],
                        [
                          4,
                          1,
                          1
                        ],
                        [
                          4,
                          2,
                          3
                        ],
                        [
                          4,
                          3,
                          1
                        ],
                        [
                          5,
                          0,
                          130
                        ],
                        [
                          5,
                          1,
                          229
                        ],
                        [
                          5,
                          2,
                          122
                        ],
                        [
                          5,
                          3,
                          69
                        ],
                        [
                          6,
                          2,
                          1
                        ],
                        [
                          6,
                          3,
                          2
                        ],
                        [
                          7,
                          0,
                          52
                        ],
                        [
                          7,
                          1,
                          80
                        ],
                        [
                          7,
                          2,
                          5
                        ],
                        [
                          7,
                          3,
                          2
                        ],
                        [
                          8,
                          1,
                          2
                        ],
                        [
                          9,
                          0,
                          3
                        ],
                        [
                          9,
                          1,
                          7
                        ],
                        [
                          9,
                          2,
                          4
                        ],
                        [
                          9,
                          3,
                          2
                        ],
                        [
                          10,
                          1,
                          1
                        ],
                        [
                          10,
                          3,
                          1
                        ],
                        [
                          11,
                          0,
                          6
                        ],
                        [
                          11,
                          1,
                          9
                        ],
                        [
                          11,
                          2,
                          4
                        ],
                        [
                          11,
                          3,
                          5
                        ],
                        [
                          12,
                          1,
                          1
                        ],
                        [
                          12,
                          2,
                          1
                        ],
                        [
                          12,
                          3,
                          2
                        ],
                        [
                          13,
                          0,
                          1
                        ],
                        [
                          13,
                          2,
                          1
                        ],
                        [
                          14,
                          1,
                          2
                        ],
                        [
                          15,
                          0,
                          1
                        ],
                        [
                          15,
                          3,
                          3
                        ],
                        [
                          16,
                          3,
                          2
                        ],
                        [
                          17,
                          1,
                          4
                        ],
                        [
                          18,
                          0,
                          1
                        ],
                        [
                          18,
                          3,
                          1
                        ],
                        [
                          19,
                          0,
                          1
                        ],
                        [
                          19,
                          1,
                          1
                        ],
                        [
                          19,
                          3,
                          1
                        ],
                        [
                          20,
                          0,
                          5
                        ],
                        [
                          20,
                          1,
                          13
                        ],
                        [
                          21,
                          0,
                          2
                        ],
                        [
                          21,
                          1,
                          3
                        ],
                        [
                          21,
                          2,
                          2
                        ],
                        [
                          21,
                          3,
                          1
                        ],
                        [
                          22,
                          0,
                          1
                        ],
                        [
                          22,
                          1,
                          2
                        ],
                        [
                          23,
                          0,
                          2
                        ],
                        [
                          23,
                          1,
                          2
                        ],
                        [
                          23,
                          2,
                          2
                        ],
                        [
                          23,
                          3,
                          1
                        ],
                        [
                          24,
                          0,
                          1
                        ],
                        [
                          24,
                          1,
                          1
                        ],
                        [
                          25,
                          1,
                          2
                        ],
                        [
                          25,
                          3,
                          1
                        ],
                        [
                          26,
                          0,
                          17
                        ],
                        [
                          26,
                          1,
                          18
                        ],
                        [
                          26,
                          2,
                          69
                        ],
                        [
                          26,
                          3,
                          64
                        ],
                        [
                          27,
                          1,
                          1
                        ],
                        [
                          27,
                          3,
                          2
                        ],
                        [
                          28,
                          0,
                          20
                        ],
                        [
                          28,
                          1,
                          29
                        ],
                        [
                          28,
                          2,
                          133
                        ],
                        [
                          28,
                          3,
                          104
                        ],
                        [
                          29,
                          0,
                          2
                        ],
                        [
                          29,
                          1,
                          5
                        ],
                        [
                          29,
                          2,
                          2
                        ],
                        [
                          29,
                          3,
                          3
                        ],
                        [
                          30,
                          0,
                          31
                        ],
                        [
                          30,
                          1,
                          48
                        ],
                        [
                          30,
                          2,
                          10
                        ],
                        [
                          30,
                          3,
                          15
                        ],
                        [
                          31,
                          0,
                          1
                        ],
                        [
                          31,
                          1,
                          2
                        ],
                        [
                          31,
                          2,
                          15
                        ],
                        [
                          31,
                          3,
                          12
                        ],
                        [
                          32,
                          0,
                          1
                        ],
                        [
                          32,
                          1,
                          1
                        ],
                        [
                          33,
                          0,
                          94
                        ],
                        [
                          33,
                          1,
                          150
                        ],
                        [
                          33,
                          2,
                          63
                        ],
                        [
                          33,
                          3,
                          39
                        ],
                        [
                          34,
                          1,
                          1
                        ],
                        [
                          34,
                          2,
                          1
                        ],
                        [
                          35,
                          1,
                          4
                        ],
                        [
                          36,
                          0,
                          1
                        ],
                        [
                          36,
                          1,
                          1
                        ],
                        [
                          37,
                          1,
                          1
                        ],
                        [
                          37,
                          3,
                          1
                        ],
                        [
                          38,
                          0,
                          1
                        ],
                        [
                          38,
                          1,
                          1
                        ],
                        [
                          39,
                          0,
                          22
                        ],
                        [
                          39,
                          1,
                          44
                        ],
                        [
                          39,
                          2,
                          1
                        ],
                        [
                          40,
                          0,
                          4
                        ],
                        [
                          40,
                          1,
                          7
                        ],
                        [
                          41,
                          0,
                          1
                        ],
                        [
                          41,
                          1,
                          2
                        ],
                        [
                          41,
                          2,
                          3
                        ],
                        [
                          42,
                          0,
                          198
                        ],
                        [
                          42,
                          1,
                          374
                        ],
                        [
                          42,
                          2,
                          181
                        ],
                        [
                          42,
                          3,
                          167
                        ],
                        [
                          43,
                          0,
                          192
                        ],
                        [
                          43,
                          1,
                          338
                        ],
                        [
                          43,
                          2,
                          5
                        ],
                        [
                          43,
                          3,
                          17
                        ],
                        [
                          44,
                          0,
                          1
                        ],
                        [
                          44,
                          1,
                          1
                        ],
                        [
                          45,
                          0,
                          1
                        ],
                        [
                          45,
                          1,
                          1
                        ],
                        [
                          45,
                          3,
                          1
                        ],
                        [
                          46,
                          0,
                          1
                        ],
                        [
                          46,
                          1,
                          1
                        ],
                        [
                          46,
                          3,
                          4
                        ],
                        [
                          47,
                          0,
                          2
                        ],
                        [
                          47,
                          1,
                          3
                        ],
                        [
                          47,
                          2,
                          1
                        ],
                        [
                          47,
                          3,
                          3
                        ],
                        [
                          48,
                          1,
                          1
                        ],
                        [
                          48,
                          2,
                          1
                        ],
                        [
                          49,
                          0,
                          2
                        ],
                        [
                          49,
                          1,
                          1
                        ],
                        [
                          50,
                          0,
                          14
                        ],
                        [
                          50,
                          1,
                          19
                        ],
                        [
                          50,
                          2,
                          6
                        ],
                        [
                          50,
                          3,
                          8
                        ],
                        [
                          51,
                          0,
                          27
                        ],
                        [
                          51,
                          1,
                          55
                        ],
                        [
                          51,
                          2,
                          1
                        ],
                        [
                          52,
                          1,
                          1
                        ],
                        [
                          52,
                          2,
                          1
                        ],
                        [
                          53,
                          2,
                          2
                        ],
                        [
                          54,
                          0,
                          9
                        ],
                        [
                          54,
                          1,
                          27
                        ],
                        [
                          54,
                          2,
                          14
                        ],
                        [
                          54,
                          3,
                          11
                        ],
                        [
                          55,
                          1,
                          1
                        ],
                        [
                          55,
                          3,
                          1
                        ],
                        [
                          56,
                          0,
                          8
                        ],
                        [
                          56,
                          1,
                          9
                        ],
                        [
                          56,
                          2,
                          2
                        ],
                        [
                          56,
                          3,
                          4
                        ],
                        [
                          57,
                          0,
                          1
                        ],
                        [
                          57,
                          1,
                          1
                        ],
                        [
                          57,
                          2,
                          1
                        ],
                        [
                          57,
                          3,
                          1
                        ],
                        [
                          58,
                          0,
                          3
                        ],
                        [
                          58,
                          1,
                          1
                        ],
                        [
                          58,
                          2,
                          1
                        ],
                        [
                          58,
                          3,
                          1
                        ],
                        [
                          59,
                          1,
                          2
                        ],
                        [
                          60,
                          0,
                          3
                        ],
                        [
                          60,
                          2,
                          1
                        ],
                        [
                          61,
                          0,
                          91
                        ],
                        [
                          61,
                          1,
                          160
                        ],
                        [
                          61,
                          2,
                          4
                        ],
                        [
                          61,
                          3,
                          3
                        ],
                        [
                          62,
                          0,
                          1
                        ],
                        [
                          62,
                          1,
                          1
                        ],
                        [
                          62,
                          2,
                          1
                        ],
                        [
                          62,
                          3,
                          2
                        ],
                        [
                          63,
                          0,
                          3
                        ],
                        [
                          63,
                          1,
                          1
                        ],
                        [
                          64,
                          0,
                          1
                        ],
                        [
                          64,
                          1,
                          1
                        ],
                        [
                          64,
                          2,
                          2
                        ],
                        [
                          64,
                          3,
                          1
                        ],
                        [
                          65,
                          2,
                          1
                        ],
                        [
                          65,
                          3,
                          1
                        ],
                        [
                          66,
                          1,
                          2
                        ],
                        [
                          66,
                          2,
                          2
                        ],
                        [
                          66,
                          3,
                          2
                        ],
                        [
                          67,
                          2,
                          1
                        ],
                        [
                          67,
                          3,
                          1
                        ],
                        [
                          68,
                          0,
                          1
                        ],
                        [
                          68,
                          1,
                          2
                        ],
                        [
                          69,
                          0,
                          1
                        ],
                        [
                          69,
                          1,
                          1
                        ]
                      ],
                      "rows": [
                        {
                          "id": "269901",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Proteobacteria",
                              "c__Gammaproteobacteria",
                              "o__Pseudomonadales",
                              "f__Pseudomonadaceae"
                            ]
                          }
                        },
                        {
                          "id": "4130483",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "137056",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "1995363",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli",
                              "o__Bacillales",
                              "f__Staphylococcaceae",
                              "g__Staphylococcus"
                            ]
                          }
                        },
                        {
                          "id": "939252",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli",
                              "o__Bacillales",
                              "f__Staphylococcaceae",
                              "g__Staphylococcus"
                            ]
                          }
                        },
                        {
                          "id": "4380971",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Clostridia",
                              "o__Clostridiales",
                              "f__Clostridiaceae",
                              "g__Clostridium",
                              "s__butyricum"
                            ]
                          }
                        },
                        {
                          "id": "1081058",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli",
                              "o__Bacillales",
                              "f__Staphylococcaceae",
                              "g__Staphylococcus"
                            ]
                          }
                        },
                        {
                          "id": "4440404",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Proteobacteria",
                              "c__Betaproteobacteria",
                              "o__Neisseriales",
                              "f__Neisseriaceae",
                              "g__Neisseria"
                            ]
                          }
                        },
                        {
                          "id": "984924",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "953758",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli",
                              "o__Bacillales",
                              "f__Staphylococcaceae",
                              "g__Staphylococcus"
                            ]
                          }
                        },
                        {
                          "id": "4316928",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "152001",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli",
                              "o__Bacillales"
                            ]
                          }
                        },
                        {
                          "id": "227083",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "4445673",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Clostridia",
                              "o__Clostridiales",
                              "f__Clostridiaceae",
                              "g__Clostridium",
                              "s__perfringens"
                            ]
                          }
                        },
                        {
                          "id": "138389",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "4427114",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "153046",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli",
                              "o__Bacillales",
                              "f__Staphylococcaceae",
                              "g__Staphylococcus"
                            ]
                          }
                        },
                        {
                          "id": "1059655",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli",
                              "o__Lactobacillales",
                              "f__Streptococcaceae",
                              "g__Streptococcus"
                            ]
                          }
                        },
                        {
                          "id": "1550056",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "979261",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "12574",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Actinobacteria",
                              "c__Actinobacteria",
                              "o__Actinomycetales",
                              "f__Actinomycetaceae",
                              "g__Actinomyces",
                              "s__"
                            ]
                          }
                        },
                        {
                          "id": "368134",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli",
                              "o__Bacillales",
                              "f__Staphylococcaceae"
                            ]
                          }
                        },
                        {
                          "id": "1039016",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "996487",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli",
                              "o__Bacillales",
                              "f__Staphylococcaceae",
                              "g__Staphylococcus"
                            ]
                          }
                        },
                        {
                          "id": "1069592",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "1112200",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "4297222",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes"
                            ]
                          }
                        },
                        {
                          "id": "923151",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "532163",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Proteobacteria",
                              "c__Alphaproteobacteria",
                              "o__Rhodobacterales",
                              "f__Rhodobacteraceae"
                            ]
                          }
                        },
                        {
                          "id": "928538",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli",
                              "o__Bacillales",
                              "f__Staphylococcaceae",
                              "g__Staphylococcus"
                            ]
                          }
                        },
                        {
                          "id": "1891556",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "114510",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Proteobacteria",
                              "c__Gammaproteobacteria",
                              "o__Enterobacteriales",
                              "f__Enterobacteriaceae"
                            ]
                          }
                        },
                        {
                          "id": "158047",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "242070",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Proteobacteria",
                              "c__Gammaproteobacteria",
                              "o__Pseudomonadales",
                              "f__Pseudomonadaceae"
                            ]
                          }
                        },
                        {
                          "id": "149265",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "919490",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "164413",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "767863",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "113773",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli",
                              "o__Lactobacillales"
                            ]
                          }
                        },
                        {
                          "id": "128604",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__[Thermi]",
                              "c__Deinococci",
                              "o__Deinococcales",
                              "f__Deinococcaceae",
                              "g__Deinococcus",
                              "s__"
                            ]
                          }
                        },
                        {
                          "id": "99882",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli",
                              "o__Lactobacillales"
                            ]
                          }
                        },
                        {
                          "id": "519673",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "630141",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "219151",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Proteobacteria",
                              "c__Gammaproteobacteria",
                              "o__Pseudomonadales",
                              "f__Moraxellaceae",
                              "g__Acinetobacter"
                            ]
                          }
                        },
                        {
                          "id": "977188",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "1121111",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "894774",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "441155",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli",
                              "o__Bacillales",
                              "f__Staphylococcaceae",
                              "g__Staphylococcus"
                            ]
                          }
                        },
                        {
                          "id": "1059977",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli",
                              "o__Lactobacillales",
                              "f__Streptococcaceae",
                              "g__Streptococcus"
                            ]
                          }
                        },
                        {
                          "id": "552922",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Proteobacteria",
                              "c__Gammaproteobacteria",
                              "o__Pseudomonadales",
                              "f__Moraxellaceae",
                              "g__Acinetobacter"
                            ]
                          }
                        },
                        {
                          "id": "2874742",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "1756274",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli",
                              "o__Lactobacillales",
                              "f__Enterococcaceae"
                            ]
                          }
                        },
                        {
                          "id": "4315958",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria"
                            ]
                          }
                        },
                        {
                          "id": "617833",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "2896107",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli",
                              "o__Bacillales",
                              "f__Staphylococcaceae"
                            ]
                          }
                        },
                        {
                          "id": "4365141",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli",
                              "o__Lactobacillales",
                              "f__Leuconostocaceae",
                              "g__Leuconostoc"
                            ]
                          }
                        },
                        {
                          "id": "356733",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "1067519",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "1068955",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "4438739",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Actinobacteria",
                              "c__Actinobacteria",
                              "o__Actinomycetales",
                              "f__Propionibacteriaceae",
                              "g__Propionibacterium"
                            ]
                          }
                        },
                        {
                          "id": "164612",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "4416988",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "1055132",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "187233",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli",
                              "o__Lactobacillales",
                              "f__Lactobacillaceae"
                            ]
                          }
                        },
                        {
                          "id": "New.CleanUp.ReferenceOTU0",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria"
                            ]
                          }
                        },
                        {
                          "id": "New.CleanUp.ReferenceOTU2",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Proteobacteria",
                              "c__Alphaproteobacteria"
                            ]
                          }
                        },
                        {
                          "id": "New.CleanUp.ReferenceOTU10",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli",
                              "o__Bacillales",
                              "f__Staphylococcaceae"
                            ]
                          }
                        },
                        {
                          "id": "New.CleanUp.ReferenceOTU27",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        },
                        {
                          "id": "New.CleanUp.ReferenceOTU36",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria"
                            ]
                          }
                        },
                        {
                          "id": "New.CleanUp.ReferenceOTU39",
                          "metadata": {
                            "taxonomy": [
                              "k__Bacteria",
                              "p__Firmicutes",
                              "c__Bacilli"
                            ]
                          }
                        }
                      ],
                      "columns": [
                        {
                          "id": "HMPMockV1.1.Even1",
                          "metadata": null
                        },
                        {
                          "id": "HMPMockV1.1.Even2",
                          "metadata": null
                        },
                        {
                          "id": "HMPMockV1.2.Staggered2",
                          "metadata": null
                        },
                        {
                          "id": "HMPMockV1.2.Staggered1",
                          "metadata": null
                        }
                      ]
                    }"""

        _mock_result_table1 = '\n'.join(
            [',Dataset,F-measure,Method,Parameters,Taxon Accuracy Rate,'
             'Taxon Detection Rate,Precision,Recall,Reference,SampleID,'
             'Spearman p,Spearman r,best N alignments,confidence,coverage,'
             'e value,e-value,max accepts,min consensus fraction,similarity',
             '402,F2,1,rdp,0,0.507167459,-0.304176543,1,1,unite-97-rep-set,'
             'm3,0.741153822,0.15430335,,0,,,,,,',
             '404,F2,1,rdp,0.1,0.507167459,-0.304176543,1,1,unite-97-rep-set,'
             'm3,0.741153822,0.15430335,,0.1,,,,,,',
             '405,F2,1,rdp,0.1,0.507167459,-0.304176543,1,1,unite-97-rep-set,'
             'm2,0.741153822,0.15430335,,0.1,,,,,,',
             '408,F2,0.933333333,rdp,0.2,0.380361021,-0.360486997,0.875,1,'
             'unite-97-rep-set,m3,0.634447233,-0.20025047,,0.2,,,,,,',
             '411,F2,0.933333333,rdp,0.3,0.380361021,-0.360486997,0.875,1,'
             'unite-97-rep-set,m3,0.634447233,-0.20025047,,0.3,,,,,,',
             '414,F2,0.933333333,rdp,0.4,0.380361021,-0.360486997,0.875,1,'
             'unite-97-rep-set,m3,0.634447233,-0.20025047,,0.4,,,,,,',
             '417,F2,0.933333333,rdp,0.5,0.380361021,-0.360486997,0.875,1,'
             'unite-97-rep-set,m3,0.634447233,-0.20025047,,0.5,,,,,,',
             '2568,F2,0.933333333,uclust,0.51:0.9:3,0.968373131,0.016870865,'
             '0.875,1,unite-97-rep-set,m3,0.394741568,0.350438322,,,,,,3,0.51,'
             '0.9',
             '2559,F2,0.875,uclust,0.51:0.8:3,0.775094897,-0.111550294,'
             '0.777777778,1,unite-97-rep-set,m3,0.705299666,0.147297986,,,,,,'
             '3,0.51,0.8',
             '2560,F2,0.875,uclust,0.51:0.8:3,0.775094897,-0.111550294,'
             '0.777777778,1,unite-97-rep-set,m2,0.705299666,0.147297986,,,,,,'
             '3,0.51,0.8',
             '6,B1,0.666666667,rdp,0.2,0.01020742,0.535556024,0.523809524,'
             '0.916666667,gg_13_8_otus,m1,0.005044461,0.57579308,,0.2,,,,,,',
             '4,B1,0.647058824,rdp,0,0.006884692,0.54720917,0.5,0.916666667,'
             'gg_13_8_otus,m1,0.003014749,0.590459586,,0,,,,,,',
             '7,B1,0.647058824,rdp,0.3,0.007585897,0.541731947,0.5,0.916666667'
             ',gg_13_8_otus,m1,0.003647366,0.580986972,,0.3,,,,,,',
             '5,B1,0.628571429,rdp,0.1,0.004629362,0.55772667,0.47826087,'
             '0.916666667,gg_13_8_otus,m1,0.001859651,0.601902881,,0.1,,,,,,',
             '8,B1,0.628571429,rdp,0.4,0.005449953,0.549139489,0.47826087,'
             '0.916666667,gg_13_8_otus,m1,0.004205291,0.562678799,,0.4,,,,,,',
             '73,B1,0.628571429,uclust,0.51:0.8:3,0.008500543,0.524553724,'
             '0.47826087,0.916666667,gg_13_8_otus,m1,0.002372594,0.590696001'
             ',,,,,,3,0.51,0.8',
             '9,B1,0.611111111,rdp,0.5,0.00572351,0.536241977,0.458333333,'
             '0.916666667,gg_13_8_otus,m1,0.002837231,0.571587285,,0.5,,,,,,',
             '76,B1,0.594594595,uclust,0.51:0.9:3,0.006571463,0.519167888,'
             '0.44,0.916666667,gg_13_8_otus,m1,0.00097195,0.608520351,,,,,,'
             '3,0.51,0.9'])

        _exp_taxa = '\n'.join(
            ['o1\tother',
             'o3\tAa;Ii;Jj;Kk;Ll;Mm;Nn',
             'o2\tAa;Bb;Cc;Dd;Ee;Ff;Hh',
             'o4\tAa;Ii;Jj;Kk;Ll;Mm;Oo'])

        _obs_taxa = '\n'.join(
            ['o2\tAa;Bb;Cc;Dd;Ee;Ff;Ii',
             'o3\tAa;Ii;Pp;Qq;Rr;Tt',
             'o4\tAa;Ii;Jj;Kk;Ll;Mm;Oo',
             'o1\tAa;Bb;Cc;Dd;Ee;Ff;Gg'])

        self.table1 = Table.from_json(json.loads(_table1))
        self.table2 = Table.from_json(json.loads(_table2))
        self.table3 = Table.from_json(json.loads(_table3))

        self.mock_result_table1 = pd.DataFrame.from_csv(
            StringIO(_mock_result_table1))

        self.tmpdir = mkdtemp()

        with open(join(self.tmpdir, 'trueish-taxonomies.tsv'), 'w') as out:
            out.write(_exp_taxa)
        with open(join(self.tmpdir, 'taxonomy.tsv'), 'w') as out:
            out.write(_obs_taxa)

    @classmethod
    def tearDownClass(self):
        rmtree(self.tmpdir)


if __name__ == "__main__":
    main()
