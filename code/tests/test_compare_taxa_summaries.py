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
from cogent.util.unit_test import TestCase, main
from taxcompare.compare_taxa_summaries import _sort_and_fill_taxa_summaries

class CompareTaxaSummariesTests(TestCase):
    """Tests for the compare_taxa_summaries.py module."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.taxa_summary1_data = (['Even1','Even2','Even3'],
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

        self.taxa_summary2_data = (['Even4','Even5','Even6'],
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

        self.taxa_summary3_data = (['Even7','Even8'], ['Eukarya'],
                                   array([[1.0, 1.0]]))

        self.taxa_summary4_data = (['Even1','Even2'], ['Eukarya'],
                                   array([[0.5, 0.6]]))

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

        obs = _sort_and_fill_taxa_summaries([self.taxa_summary1_data,
                                             self.taxa_summary2_data,
                                             self.taxa_summary3_data])
        self.assertFloatEqual(obs, exp)

    def test_sort_and_fill_taxa_summaries_same(self):
        """Test when the taxa summaries are the same."""
        exp = [(['Even7','Even8'], ['Eukarya'], array([[1.0, 1.0]])),
               (['Even7','Even8'], ['Eukarya'], array([[1.0, 1.0]]))]
        obs = _sort_and_fill_taxa_summaries([self.taxa_summary3_data,
                                            self.taxa_summary3_data])
        self.assertFloatEqual(obs, exp)

        exp = [(['Even7','Even8'], ['Eukarya'], array([[1.0, 1.0]])),
               (['Even1','Even2'], ['Eukarya'], array([[0.5, 0.6]]))]
        obs = _sort_and_fill_taxa_summaries([self.taxa_summary3_data,
                                            self.taxa_summary4_data])
        self.assertFloatEqual(obs, exp)


if __name__ == "__main__":
    main()
