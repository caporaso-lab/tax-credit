#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout", "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"
__status__ = "Development"

"""Test suite for the compare_taxa_summaries.py module."""

from numpy import array
from cogent.maths.stats.test import pearson
from cogent.util.unit_test import TestCase, main
from taxcompare.compare_taxa_summaries import (_compare_all_to_expected,
        compare_taxa_summaries, _sort_and_fill_taxa_summaries,
        spearman_correlation)

class CompareTaxaSummariesTests(TestCase):
    """Tests for the compare_taxa_summaries.py module."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.taxa_summary1 = (['Even1','Even2','Even3'],
            ['Bacteria;Actinobacteria;Actinobacteria(class);Actinobacteridae',
             'Bacteria;Bacteroidetes/Chlorobigroup;Bacteroidetes;Bacteroidia',
             'Bacteria;Firmicutes;Bacilli;Lactobacillales',
             'Bacteria;Firmicutes;Clostridia;Clostridiales',
             'Bacteria;Firmicutes;Erysipelotrichi;Erysipelotrichales',
             'Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales',
             'No blast hit;Other'],
            array([[0.0880247251673, 0.0721968465746, 0.081371761759],
                   [0.192137761955, 0.191095101593, 0.188504131885],
                   [0.0264895739603, 0.0259942669171, 0.0318460745596],
                   [0.491800007824, 0.526186212556, 0.49911159984],
                   [0.0311411916592, 0.0184083913576, 0.0282325481054],
                   [0.166137214246, 0.163087129528, 0.168923372865],
                   [0.00426952518811, 0.00303205147361, 0.0020105109874]]))

        self.taxa_summary2 = (['Even4','Even5','Even6'],
            ['Bacteria;Actinobacteria;Actinobacteria(class);NotARealTaxa',
             'Bacteria;AnotherFakeTaxa',
             'Bacteria;Bacteroidetes/Chlorobigroup;Bacteroidetes;Bacteroidia',
             'Bacteria;Firmicutes;Bacilli;Lactobacillales',
             'Bacteria;Firmicutes;Clostridia;Clostridiales',
             'Bacteria;Firmicutes;Erysipelotrichi;Erysipelotrichales',
             'Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales',
             'No blast hit;Other'],
            array([[0.99, 0.11, 0.075],
                   [0.1921, 0.19109, 0.18],
                   [0.192137761955, 0.191095101593, 0.188504131885],
                   [0.0264895739603, 0.0259942669171, 0.0318460745596],
                   [0.491800007824, 0.526186212556, 0.49911159984],
                   [0.0311411916592, 0.0184083913576, 0.0282325481054],
                   [0.166137214246, 0.163087129528, 0.168923372865],
                   [0.00426952518811, 0.00303205147361, 0.0020105109874]]))

        self.taxa_summary3 = (['Even7','Even8'], ['Eukarya'],
                                   array([[1.0, 1.0]]))

        # Has the same taxa as self.taxa_summary3_data (above).
        self.taxa_summary4 = (['Even1','Even2'], ['Eukarya'],
                                   array([[0.5, 0.6]]))

        # No intersection with self.taxa_summary4_data.
        self.taxa_summary5 = (['Even1','Even2'], ['foo'],
                                   array([[0.5, 0.6]]))

        # For testing all versus expected comparison mode.
        self.taxa_summary_exp1 = (['Expected'], ['Eukarya', 'Bacteria'],
                                   array([[0.5], [0.6]]))
        self.taxa_summary_obs1 = (['S1', 'S2'], ['Eukarya', 'Bacteria'],
                                   array([[0.4, 0.5], [0.5, 0.7]]))
        # For testing all versus expected using spearman (contains repeats).
        self.taxa_summary_exp2 = (['Expected'],
                                  ['Eukarya', 'Bacteria', 'Archaea'],
                                   array([[0.5], [0.6], [0.4]]))
        self.taxa_summary_obs2 = (['S1', 'S2'], ['Eukarya', 'Bacteria'],
                                   array([[0.4, 0.5], [0.5, 0.7], [0.4, 0.4]]))

    def test_compare_taxa_summaries_invalid_comparison_mode(self):
        """Test with an invalid (unrecognized) comparison mode."""
        self.assertRaises(ValueError, compare_taxa_summaries,
                self.taxa_summary1, self.taxa_summary5, 'cool comparison')

    def test_compare_taxa_summaries_invalid_correlation_type(self):
        """Test with an invalid (unrecognized) correlation type."""
        self.assertRaises(ValueError, compare_taxa_summaries,
                self.taxa_summary1, self.taxa_summary5, 'paired', 'foo')

    def test_sort_and_fill_taxa_summaries(self):
        """Test _sort_and_fill_taxa_summaries functions as expected."""
        exp = [
            (['Even1','Even2','Even3'],
             ['Bacteria;Actinobacteria;Actinobacteria(class);Actinobacteridae',
              'Bacteria;Actinobacteria;Actinobacteria(class);NotARealTaxa',
              'Bacteria;AnotherFakeTaxa',
              'Bacteria;Bacteroidetes/Chlorobigroup;Bacteroidetes;Bacteroidia',
              'Bacteria;Firmicutes;Bacilli;Lactobacillales',
              'Bacteria;Firmicutes;Clostridia;Clostridiales',
              'Bacteria;Firmicutes;Erysipelotrichi;Erysipelotrichales',
              'Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales',
              'Eukarya',
              'No blast hit;Other'],
             array([[0.0880247251673, 0.0721968465746, 0.081371761759],
                    [0.,0.,0.],
                    [0.,0.,0.],
                    [0.192137761955, 0.191095101593, 0.188504131885],
                    [0.0264895739603, 0.0259942669171, 0.0318460745596],
                    [0.491800007824, 0.526186212556, 0.49911159984],
                    [0.0311411916592, 0.0184083913576, 0.0282325481054],
                    [0.166137214246, 0.163087129528, 0.168923372865],
                    [0.,0.,0.],
                    [0.00426952518811, 0.00303205147361, 0.0020105109874]])),
            (['Even4','Even5','Even6'],
             ['Bacteria;Actinobacteria;Actinobacteria(class);Actinobacteridae',
              'Bacteria;Actinobacteria;Actinobacteria(class);NotARealTaxa',
              'Bacteria;AnotherFakeTaxa',
              'Bacteria;Bacteroidetes/Chlorobigroup;Bacteroidetes;Bacteroidia',
              'Bacteria;Firmicutes;Bacilli;Lactobacillales',
              'Bacteria;Firmicutes;Clostridia;Clostridiales',
              'Bacteria;Firmicutes;Erysipelotrichi;Erysipelotrichales',
              'Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales',
              'Eukarya',
              'No blast hit;Other'],
             array([[0.,0.,0.],
                    [0.99, 0.11, 0.075],
                    [0.1921, 0.19109, 0.18],
                    [0.192137761955, 0.191095101593, 0.188504131885],
                    [0.0264895739603, 0.0259942669171, 0.0318460745596],
                    [0.491800007824, 0.526186212556, 0.49911159984],
                    [0.0311411916592, 0.0184083913576, 0.0282325481054],
                    [0.166137214246, 0.163087129528, 0.168923372865],
                    [0.,0.,0.],
                    [0.00426952518811, 0.00303205147361, 0.0020105109874]])),
            (['Even7','Even8'],
             ['Bacteria;Actinobacteria;Actinobacteria(class);Actinobacteridae',
              'Bacteria;Actinobacteria;Actinobacteria(class);NotARealTaxa',
              'Bacteria;AnotherFakeTaxa',
              'Bacteria;Bacteroidetes/Chlorobigroup;Bacteroidetes;Bacteroidia',
              'Bacteria;Firmicutes;Bacilli;Lactobacillales',
              'Bacteria;Firmicutes;Clostridia;Clostridiales',
              'Bacteria;Firmicutes;Erysipelotrichi;Erysipelotrichales',
              'Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales',
              'Eukarya',
              'No blast hit;Other'],
             array([[0.,0.],
                    [0.,0.],
                    [0.,0.],
                    [0.,0.],
                    [0.,0.],
                    [0.,0.],
                    [0.,0.],
                    [0.,0.],
                    [1.,1.],
                    [0.,0.]]))
        ]

        obs = _sort_and_fill_taxa_summaries([self.taxa_summary1,
                                             self.taxa_summary2,
                                             self.taxa_summary3])
        self.assertFloatEqual(obs, exp)

    def test_sort_and_fill_taxa_summaries_same(self):
        """Test when the taxa summaries are the same."""
        exp = [(['Even7','Even8'], ['Eukarya'], array([[1.0, 1.0]])),
               (['Even7','Even8'], ['Eukarya'], array([[1.0, 1.0]]))]
        obs = _sort_and_fill_taxa_summaries([self.taxa_summary3,
                                             self.taxa_summary3])
        self.assertFloatEqual(obs, exp)

        exp = [(['Even7','Even8'], ['Eukarya'], array([[1.0, 1.0]])),
               (['Even1','Even2'], ['Eukarya'], array([[0.5, 0.6]]))]
        obs = _sort_and_fill_taxa_summaries([self.taxa_summary3,
                                             self.taxa_summary4])
        self.assertFloatEqual(obs, exp)

        # Test the other direction.
        exp = [(['Even1','Even2'], ['Eukarya'], array([[0.5, 0.6]])),
               (['Even7','Even8'], ['Eukarya'], array([[1.0, 1.0]]))]
        obs = _sort_and_fill_taxa_summaries([self.taxa_summary4,
                                             self.taxa_summary3])
        self.assertFloatEqual(obs, exp)

    def test_sort_and_fill_taxa_summaries_no_intersection(self):
        """Test when the taxa summaries have no intersection."""
        exp = [(['Even1','Even2'], ['Eukarya', 'foo'], array([[0.5, 0.6],
                                                              [0.0, 0.0]])),
               (['Even1','Even2'], ['Eukarya', 'foo'], array([[0.0, 0.0],
                                                              [0.5, 0.6]]))]
        obs = _sort_and_fill_taxa_summaries([self.taxa_summary4,
                                             self.taxa_summary5])
        self.assertFloatEqual(obs, exp)

    def test_compare_all_to_expected_pearson(self):
        """Test _compare_all_to_expected() functions correctly with pearson."""
        exp = [1.0, 1.0]
        obs = _compare_all_to_expected(self.taxa_summary_obs1,
                                       self.taxa_summary_exp1, pearson)
        self.assertFloatEqual(obs, exp)

    def test_compare_all_to_expected_spearman(self):
        """Test _compare_all_to_expected() functions with spearman."""
        # No repeats in ranks.
        exp = [1.0, 1.0]
        obs = _compare_all_to_expected(self.taxa_summary_obs1,
                                       self.taxa_summary_exp1,
                                       spearman_correlation)
        self.assertFloatEqual(obs, exp)

        # Repeats in ranks.
        exp = [0.866025, 1.0]
        obs = _compare_all_to_expected(self.taxa_summary_obs2,
                                       self.taxa_summary_exp2,
                                       spearman_correlation)
        self.assertFloatEqual(obs, exp)

    def test_compare_all_to_expected_invalid_expected_sample_count(self):
        """Test running on expected taxa summary without exactly one sample."""
        self.assertRaises(ValueError, _compare_all_to_expected,
                          self.taxa_summary1, self.taxa_summary2, pearson)


if __name__ == "__main__":
    main()
