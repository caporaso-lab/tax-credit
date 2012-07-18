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

"""Test suite for the compare_all_to_expected.py module."""

from numpy import array
from cogent.util.unit_test import TestCase, main
from taxcompare.compare_taxa_summaries import (add_filename_suffix,
        compare_all_to_expected, _compute_all_to_expected_correlations,
        format_taxa_summary, _get_correlation_function,
        _make_compatible_taxa_summaries, _sort_and_fill_taxa_summaries,
        _pearson_correlation, _spearman_correlation)

class CompareTaxaSummariesTests(TestCase):
    """Tests for the compare_all_to_expected.py module."""

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

        # Different sample ID from self.taxa_summary5.
        self.taxa_summary6 = (['Even1','Even3'], ['foo'],
                                   array([[0.1, 0.6]]))

        # Samples are in different orders from self.taxa_summary6.
        self.taxa_summary7 = (['Even3','Even1'], ['foo'],
                                   array([[0.2, 0.77]]))

        # Samples are not in alphabetical order, and we have multiple taxa.
        self.taxa_summary8 = (['Even3','Even1', 'S7'], ['foo', 'bar'],
                                   array([[0.2, 0.77, 0.001],
                                          [0.45, 0.9, 0.0]]))

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

    def test_compare_all_to_expected_pearson(self):
        """Functions correctly using 'expected' comparison mode and pearson."""
        exp = ((['S1', 'S2'], ['Bacteria', 'Eukarya'], array([[ 0.5,  0.7],
            [ 0.4,  0.5]])), (['Expected'], ['Bacteria', 'Eukarya'],
            array([[ 0.6], [ 0.5]])), '# All samples in the first taxa '
            'summary file are compared to the single sample "Expected" in the '
            'second taxa summary file using pearson correlation.\n\tExpected\n'
            'S1\t1.0\nS2\t1.0')

        obs = compare_all_to_expected(self.taxa_summary_obs1,
                self.taxa_summary_exp1, 'pearson')
        self.assertEqual(obs, exp)

    def test_compare_all_to_expected_expected_spearman(self):
        """Functions correctly using 'expected' comparison mode w/ spearman."""
        exp = ((['S1', 'S2'], ['Bacteria', 'Eukarya'], array([[ 0.5,  0.7],
            [ 0.4,  0.5]])), (['Expected'], ['Bacteria', 'Eukarya'],
            array([[ 0.6], [ 0.5]])), '# All samples in the first taxa '
            'summary file are compared to the single sample "Expected" in the '
            'second taxa summary file using spearman correlation.\n\tExpected'
            '\nS1\t1.0\nS2\t1.0')

        obs = compare_all_to_expected(self.taxa_summary_obs1,
                self.taxa_summary_exp1, 'spearman')
        self.assertEqual(obs, exp)

    def test_add_filename_suffix(self):
        """Test adding a suffix to a filename works correctly."""
        self.assertEqual(add_filename_suffix('/foo/bar/baz.txt', 'z'),
                                             'bazz.txt')
        self.assertEqual(add_filename_suffix('baz.txt', 'z'),
                                             'bazz.txt')
        self.assertEqual(add_filename_suffix('/foo/bar/baz', 'z'),
                                             'bazz')
        self.assertEqual(add_filename_suffix('baz', 'z'),
                                             'bazz')
        self.assertEqual(add_filename_suffix('/baz.fasta.txt', 'z'),
                                             'baz.fastaz.txt')
        self.assertEqual(add_filename_suffix('baz.fasta.txt', 'z'),
                                             'baz.fastaz.txt')
        self.assertEqual(add_filename_suffix('/foo/', 'z'), 'z')

    def test_format_taxa_summary(self):
        """Test formatting a taxa summary works correctly."""
        # More than one sample.
        exp = 'Taxon\tEven7\tEven8\nEukarya\t1.0\t1.0\n'
        obs = format_taxa_summary(self.taxa_summary3)
        self.assertEqual(obs, exp)

        # More than one taxon.
        exp = 'Taxon\tExpected\nEukarya\t0.5\nBacteria\t0.6\nArchaea\t0.4\n'
        obs = format_taxa_summary(self.taxa_summary_exp2)
        self.assertEqual(obs, exp)

    def test_get_correlation_function(self):
        """Test returns correct correlation function."""
        self.assertEqual(_get_correlation_function('pearson'),
                         _pearson_correlation)
        self.assertEqual(_get_correlation_function('spearman'),
                         _spearman_correlation)

    def test_get_correlation_function_invalid_correlation_type(self):
        """Test with an invalid (unrecognized) correlation type."""
        self.assertRaises(ValueError, _get_correlation_function, 'foo')

    def test_make_compatible_taxa_summaries_one_ts(self):
        """Test making compatible taxa summaries works correctly on one ts."""
        # Should return the same taxa summary since no samples should get
        # stripped out, and the sample IDs are already in alphabetical order.
        obs = _make_compatible_taxa_summaries([self.taxa_summary1])
        self.assertFloatEqual(obs, [self.taxa_summary1])

        obs = _make_compatible_taxa_summaries([self.taxa_summary6])
        self.assertFloatEqual(obs, [self.taxa_summary6])

        # The sample IDs will be reordered to be in alphabetical order.
        exp = [(['Even1', 'Even3'], ['foo'], array([[ 0.77,  0.2 ]]))]
        obs = _make_compatible_taxa_summaries([self.taxa_summary7])
        self.assertFloatEqual(obs, exp)

        # The sample IDs will be reordered to be in alphabetical order (this
        # one includes multiple taxa).
        exp = [(['Even1', 'Even3', 'S7'], ['foo', 'bar'],
                array([[ 0.77 ,  0.2  ,  0.001], [ 0.9  ,  0.45 ,  0.   ]]))]
        obs = _make_compatible_taxa_summaries([self.taxa_summary8])
        self.assertFloatEqual(obs, exp)

    def test_make_compatible_taxa_summaries_two_ts(self):
        """Test making compatible taxa summaries works correctly on two ts."""
        exp = [(['Even1'], ['foo'], array([[ 0.5]])), (['Even1'], ['foo'],
                array([[ 0.1]]))]
        obs = _make_compatible_taxa_summaries([self.taxa_summary5,
            self.taxa_summary6])
        self.assertFloatEqual(obs, exp)

    def test_make_compatible_taxa_summaries_multiple_ts(self):
        """Test making compatible taxa summaries works correctly on many ts."""
        exp = [(['Even1'], ['foo'], array([[ 0.5]])), (['Even1'], ['foo'],
                array([[ 0.1]])), (['Even1'], ['Eukarya'], array([[ 0.5]]))]
        obs = _make_compatible_taxa_summaries([self.taxa_summary5,
            self.taxa_summary6, self.taxa_summary4])
        self.assertFloatEqual(obs, exp)

    def test_make_compatible_taxa_summaries_already_compatible(self):
        """Test on taxa summaries that are already compatible."""
        # Using two compatible ts.
        obs = _make_compatible_taxa_summaries([self.taxa_summary4,
                                               self.taxa_summary5])
        self.assertFloatEqual(obs, [self.taxa_summary4, self.taxa_summary5])

        # Using four compatible ts.
        obs = _make_compatible_taxa_summaries([self.taxa_summary4,
                                               self.taxa_summary5,
                                               self.taxa_summary5,
                                               self.taxa_summary4])
        self.assertFloatEqual(obs, [self.taxa_summary4, self.taxa_summary5,
                                    self.taxa_summary5, self.taxa_summary4])

    def test_make_compatible_taxa_summaries_incompatible(self):
        """Test on taxa summaries that have no common sample IDs."""
        # Using two incompatible ts.
        self.assertRaises(ValueError, _make_compatible_taxa_summaries,
                          [self.taxa_summary3, self.taxa_summary4])

        # Using three incompatible ts.
        self.assertRaises(ValueError, _make_compatible_taxa_summaries,
            [self.taxa_summary1, self.taxa_summary2, self.taxa_summary3])

    def test_make_compatible_taxa_summaries_unordered(self):
        """Test on taxa summaries whose sample IDs are in different orders."""
        # Using two compatible ts.
        exp = [(['Even1', 'Even3'], ['foo'], array([[ 0.1,  0.6]])),
               (['Even1', 'Even3'], ['foo'], array([[ 0.77,  0.2 ]]))]
        obs = _make_compatible_taxa_summaries(
                [self.taxa_summary6, self.taxa_summary7])
        self.assertFloatEqual(obs, exp)

        #self.taxa_summary6 = (['Even1','Even3'], ['foo'],
        #                           array([[0.1, 0.6]]))
        #self.taxa_summary7 = (['Even3','Even1'], ['foo'],
        #                           array([[0.2, 0.77]]))

        # Using three incompatible ts.
        #self.assertRaises(ValueError, _make_compatible_taxa_summaries,
        #    [self.taxa_summary1, self.taxa_summary2, self.taxa_summary3])

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

    def test_compute_all_to_expected_correlations_pearson(self):
        """Test functions correctly with pearson correlation."""
        exp = [1.0, 1.0]
        obs = _compute_all_to_expected_correlations(self.taxa_summary_obs1,
                                       self.taxa_summary_exp1,
                                       _pearson_correlation)
        self.assertFloatEqual(obs, exp)

    def test_compute_all_to_expected_correlations_spearman(self):
        """Test functions correctly with spearman correlation."""
        # No repeats in ranks.
        exp = [1.0, 1.0]
        obs = _compute_all_to_expected_correlations(self.taxa_summary_obs1,
                                       self.taxa_summary_exp1,
                                       _spearman_correlation)
        self.assertFloatEqual(obs, exp)

        # Repeats in ranks.
        exp = [0.866025, 1.0]
        obs = _compute_all_to_expected_correlations(self.taxa_summary_obs2,
                                       self.taxa_summary_exp2,
                                       _spearman_correlation)
        self.assertFloatEqual(obs, exp)

    def test_compute_all_to_expected_correlations_invalid_sample_count(self):
        """Test running on expected taxa summary without exactly one sample."""
        self.assertRaises(ValueError, _compute_all_to_expected_correlations,
                          self.taxa_summary1, self.taxa_summary2,
                          _pearson_correlation)

    def test_pearson_correlation_invalid_input(self):
        """Test running pearson correlation on bad input."""
        self.assertRaises(ValueError, _pearson_correlation,
                          [1.4, 2.5], [5.6, 8.8, 9.0])
        self.assertRaises(ValueError, _pearson_correlation, [1.4], [5.6])


if __name__ == "__main__":
    main()
