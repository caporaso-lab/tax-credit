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

from taxcompare.multiple_assign_taxonomy import (
        assign_taxonomy_multiple_times, _generate_rdp_commands,
        _generate_taxa_processing_commands)

class MultipleAssignTaxonomyTests(TestCase):
    """Tests for the multiple_assign_taxonomy.py module."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        pass

    def test_assign_taxonomy_multiple_times(self):
        """Functions correctly using standard valid input data."""
        pass

    def test_generate_rdp_commands(self):
        """Functions correctly using standard valid input data."""
        exp = [[('Assigning taxonomy (RDP, 0.8 confidence)',
            'assign_taxonomy.py -i /foo/bar/rep_set.fna -o /foo/bar/rdp_0.8 '
            '-c 0.8 -m rdp -r /baz/reference_seqs.fasta -t '
            '/baz/rdp_id_to_taxonomy.txt')],
            [('Adding taxa (RDP, 0.8 confidence)',
            'add_taxa.py -i /foo/bar/otu_table.biom -o '
            '/foo/bar/rdp_0.8/otu_table_w_taxa.biom -t '
            '/foo/bar/rdp_0.8/rep_set_tax_assignments.txt')],
            [('Summarizing taxa (RDP, 0.8 confidence)',
            'summarize_taxa.py -i /foo/bar/rdp_0.8/otu_table_w_taxa.biom -o '
            '/foo/bar/rdp_0.8')],
            [('Assigning taxonomy (RDP, 0.6 confidence)',
            'assign_taxonomy.py -i /foo/bar/rep_set.fna -o /foo/bar/rdp_0.6 '
            '-c 0.6 -m rdp -r /baz/reference_seqs.fasta -t '
            '/baz/rdp_id_to_taxonomy.txt')],
            [('Adding taxa (RDP, 0.6 confidence)',
            'add_taxa.py -i /foo/bar/otu_table.biom -o '
            '/foo/bar/rdp_0.6/otu_table_w_taxa.biom -t '
            '/foo/bar/rdp_0.6/rep_set_tax_assignments.txt')],
            [('Summarizing taxa (RDP, 0.6 confidence)',
            'summarize_taxa.py -i /foo/bar/rdp_0.6/otu_table_w_taxa.biom -o '
            '/foo/bar/rdp_0.6')]]

        obs = _generate_rdp_commands('/foo/bar', '/foo/bar/rep_set.fna',
                '/baz/reference_seqs.fasta', '/baz/rdp_id_to_taxonomy.txt',
                '/foo/bar/otu_table.biom', [0.80, 0.60])
        self.assertEqual(obs, exp)

    def test_generate_taxa_processing_commands(self):
        """Functions correctly using standard valid input data."""
        exp = ([('Adding taxa (RDP, 0.8 confidence)',
            'add_taxa.py -i /foo/otu_table.biom -o '
            '/foo/rdp_0.8/otu_table_w_taxa.biom -t '
            '/foo/rdp_0.8/rep_set_tax_assignments.txt')],
            [('Summarizing taxa (RDP, 0.8 confidence)',
            'summarize_taxa.py -i /foo/rdp_0.8/otu_table_w_taxa.biom -o '
            '/foo/rdp_0.8')])
        obs = _generate_taxa_processing_commands('/foo/rdp_0.8',
                '/foo/rep_set.fna', '/foo/otu_table.biom',
                'RDP, 0.8 confidence')
        self.assertEqual(obs, exp)


if __name__ == "__main__":
    main()
