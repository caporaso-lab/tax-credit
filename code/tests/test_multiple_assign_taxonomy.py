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

from os import makedirs, getcwd, chdir
from os.path import exists
from shutil import rmtree
from tempfile import mkdtemp, NamedTemporaryFile
from cogent.util.misc import remove_files
from cogent.util.unit_test import TestCase, main
from qiime.test import initiate_timeout, disable_timeout
from qiime.util import get_qiime_temp_dir, get_tmp_filename
from qiime.workflow import WorkflowError

from taxcompare.multiple_assign_taxonomy import (
        assign_taxonomy_multiple_times,
        _generate_rdp_commands,
        _generate_blast_commands,
        _generate_mothur_commands,
        _generate_rtax_commands,
        _generate_taxa_processing_commands)

class MultipleAssignTaxonomyTests(TestCase):
    """Tests for the multiple_assign_taxonomy.py module."""

    def setUp(self):
        """Set up files/environment that will be used by the tests."""
        # The prefix to use for temporary files. This prefix may be added to,
        # but all temp dirs and files created by the tests will have this
        # prefix at a minimum.
        self.prefix = 'multiple_assign_taxonomy_tests'

        self.start_dir = getcwd()
        self.dirs_to_remove = []
        self.files_to_remove = []
        
        self.tmp_dir = get_qiime_temp_dir()
        if not exists(self.tmp_dir):
            makedirs(self.tmp_dir)
            # if test creates the temp dir, also remove it
            self.dirs_to_remove.append(self.tmp_dir)

        # setup temporary output directories
        self.output_dir = mkdtemp(dir=self.tmp_dir,
                             prefix='%s_output_dir_' %self.prefix)
        self.dirs_to_remove.append(self.output_dir)

        initiate_timeout(60)

    def tearDown(self):
        """ """
        disable_timeout()

        # change back to the start dir - some workflows change directory
        chdir(self.start_dir)

        remove_files(self.files_to_remove)
        # remove directories last, so we don't get errors
        # trying to remove files which may be in the directories
        for d in self.dirs_to_remove:
            if exists(d):
                rmtree(d)

    def test_assign_taxonomy_multiple_times(self):
        """Functions correctly using standard valid input data."""
        pass

    def test_assign_taxonomy_multiple_times_invalid_input(self):
        """Test that errors are thrown using various types of invalid input."""
        # link local output_dir to global output_dir
        out_dir = self.output_dir
        # The output directory already exists, and we aren't in force mode.
        self.assertRaises(WorkflowError, assign_taxonomy_multiple_times,
                ['/foo', '/bar/'], out_dir, ['rdp', 'blast'],
                '/foo/ref_seqs.fasta', 'in.fasta', 'otu.biom')

        # The input directories don't exist.
        self.assertRaises(WorkflowError, assign_taxonomy_multiple_times,
                ['/foobarbaz', '/foobarbaz2/'], out_dir, ['rdp', 'blast'],
                '/foo/ref_seqs.fasta', 'in.fasta', 'otu.biom', force=True)

        # Invalid assignment method.
        self.assertRaises(WorkflowError, assign_taxonomy_multiple_times,
                          [out_dir, out_dir], out_dir, ['foo', 'rdp'],
                '/foo/ref_seqs.fasta', 'in.fasta', 'otu.biom', force=True)

        # ID to taxonomy map is missing
        self.assertRaises(WorkflowError, assign_taxonomy_multiple_times,
                          [out_dir, out_dir], out_dir, ['rdp', 'mothur'],
                          '/foo/ref_seqs.fasta', 'in.fasta', 'otu.biom',
                          force=True)

    # test bad RDP input
    def test_invalid_rdp_input(self):
        """Test that errors are thrown using invalid input for RDP."""
        # link local output_dir to global output_dir
        out_dir = self.output_dir
        # RDP confidences are missing.
        self.assertRaises(WorkflowError, assign_taxonomy_multiple_times,
                          [out_dir, out_dir], out_dir, ['rdp', 'rdp'],
                          '/foo/ref_seqs.fasta', 'in.fasta', 'otu.biom',
                          id_to_taxonomy_fp='/foo/id_to_tax.txt', force=True)

    # test bad BLAST input
    def test_invalid_blast_input(self):
        """Test that errors are thrown using invalid input for BLAST."""
        # link local output_dir to global output_dir
        out_dir = self.output_dir
        # BLAST E-values are missing.
        self.assertRaises(WorkflowError, assign_taxonomy_multiple_times,
                          [out_dir, out_dir], out_dir, ['blast', 'blast'],
                          '/foo/ref_seqs.fasta', 'in.fasta', 'otu.biom',
                          id_to_taxonomy_fp='/foo/id_to_tax.txt', force=True)

    # test bad Mothur input
    def test_invalid_mothur_input(self):
        """Test that errors are thrown using invalid input for Mothur."""
        # link local output_dir to global output_dir
        out_dir = self.output_dir
        # Mothur confidences are missing.
        self.assertRaises(WorkflowError, assign_taxonomy_multiple_times,
                          [out_dir, out_dir], out_dir, ['mothur', 'mothur'],
                          '/foo/ref_seqs.fasta', 'in.fasta', 'otu.biom',
                          id_to_taxonomy_fp='/foo/id_to_tax.txt', force=True)

    # test bad RTAX input
    def test_invalid_rtax_input(self):
        """Test that errors are thrown using invalid input for RTAX."""
        # link local output_dir to global output_dir
        out_dir = self.output_dir
        # RTAX first read is missing.
        self.assertRaises(WorkflowError, assign_taxonomy_multiple_times,
                          [out_dir, out_dir], out_dir, ['rtax', 'rtax'],
                          '/foo/ref_seqs/fasta', 'in.fasta', 'otu.biom',
                          id_to_taxonomy_fp='/foo/id_to_tax.txt', force=True)

    # test RDP command generation
    def test_generate_rdp_commands(self):
        """Functions correctly using standard valid input data."""
        exp = [[('Assigning taxonomy (RDP, 0.8 confidence)',
                 'assign_taxonomy.py -i /foo/bar/rep_set.fna -o /foo/bar/rdp_0.8 '
                 '-c 0.8 -m rdp -r /baz/reference_seqs.fasta -t /baz/id_to_taxonomy.txt')],
               [('Adding taxa (RDP, 0.8 confidence)',
                 'add_taxa.py -i /foo/bar/otu_table.biom -o '
                 '/foo/bar/rdp_0.8/otu_table_w_taxa.biom -t '
                 '/foo/bar/rdp_0.8/rep_set_tax_assignments.txt')],
               [('Summarizing taxa (RDP, 0.8 confidence)',
                 'summarize_taxa.py -i /foo/bar/rdp_0.8/otu_table_w_taxa.biom -o '
                 '/foo/bar/rdp_0.8')],
               [('Assigning taxonomy (RDP, 0.6 confidence)',
                 'assign_taxonomy.py -i /foo/bar/rep_set.fna -o /foo/bar/rdp_0.6 '
                 '-c 0.6 -m rdp -r /baz/reference_seqs.fasta -t /baz/id_to_taxonomy.txt')],
               [('Adding taxa (RDP, 0.6 confidence)',
                 'add_taxa.py -i /foo/bar/otu_table.biom -o '
                 '/foo/bar/rdp_0.6/otu_table_w_taxa.biom -t '
                 '/foo/bar/rdp_0.6/rep_set_tax_assignments.txt')],
               [('Summarizing taxa (RDP, 0.6 confidence)',
                 'summarize_taxa.py -i /foo/bar/rdp_0.6/otu_table_w_taxa.biom -o '
                 '/foo/bar/rdp_0.6')]]

        obs = _generate_rdp_commands('/foo/bar', '/foo/bar/rep_set.fna',
                '/baz/reference_seqs.fasta', '/baz/id_to_taxonomy.txt',
                '/foo/bar/otu_table.biom', [0.80, 0.60])
        self.assertEqual(obs, exp)

    # test blast command generation
    def test_generate_blast_commands(self):
        """Functions correctly using standard valid input data."""
        exp = [[('Assigning taxonomy (BLAST, E 0.002)',
                 'assign_taxonomy.py -i /foo/bar/rep_set.fna -o /foo/bar/blast_0.002 -e 0.002 '
                 '-m blast -r /baz/reference_seqs.fasta -t /baz/id_to_taxonomy.txt')],
               [('Adding taxa (BLAST, E 0.002)',
                 'add_taxa.py -i /foo/bar/otu_table.biom -o /foo/bar/blast_0.002/otu_table_w_taxa.biom '
                 '-t /foo/bar/blast_0.002/rep_set_tax_assignments.txt')],
               [('Summarizing taxa (BLAST, E 0.002)',
                 'summarize_taxa.py -i /foo/bar/blast_0.002/otu_table_w_taxa.biom '
                 '-o /foo/bar/blast_0.002')],
               [('Assigning taxonomy (BLAST, E 0.005)',
                 'assign_taxonomy.py -i /foo/bar/rep_set.fna -o /foo/bar/blast_0.005 -e 0.005 '
                 '-m blast -r /baz/reference_seqs.fasta -t /baz/id_to_taxonomy.txt')],
               [('Adding taxa (BLAST, E 0.005)',
                 'add_taxa.py -i /foo/bar/otu_table.biom -o /foo/bar/blast_0.005/otu_table_w_taxa.biom '
                 '-t /foo/bar/blast_0.005/rep_set_tax_assignments.txt')],
               [('Summarizing taxa (BLAST, E 0.005)',
                 'summarize_taxa.py -i /foo/bar/blast_0.005/otu_table_w_taxa.biom '
                 '-o /foo/bar/blast_0.005')]]

        obs = _generate_blast_commands('/foo/bar', '/foo/bar/rep_set.fna',
                '/baz/reference_seqs.fasta', '/baz/id_to_taxonomy.txt',
                '/foo/bar/otu_table.biom', [0.002, 0.005])
        self.assertEqual(obs, exp)

    # test mothur command generation
    def test_generate_mothur_commands(self):
        """Functions correctly using standard valid input data."""
        exp = [[('Assigning taxonomy (Mothur, 0.8 confidence)',
                 'assign_taxonomy.py -i /foo/bar/rep_set.fna -o /foo/bar/mothur_0.8 '
                 '-c 0.8 -m mothur -r /baz/reference_seqs.fasta -t '
                 '/baz/id_to_taxonomy.txt')],
               [('Adding taxa (Mothur, 0.8 confidence)',
                 'add_taxa.py -i /foo/bar/otu_table.biom -o '
                 '/foo/bar/mothur_0.8/otu_table_w_taxa.biom -t '
                 '/foo/bar/mothur_0.8/rep_set_tax_assignments.txt')],
               [('Summarizing taxa (Mothur, 0.8 confidence)',
                 'summarize_taxa.py -i /foo/bar/mothur_0.8/otu_table_w_taxa.biom -o '
                 '/foo/bar/mothur_0.8')],
               [('Assigning taxonomy (Mothur, 0.6 confidence)',
                 'assign_taxonomy.py -i /foo/bar/rep_set.fna -o /foo/bar/mothur_0.6 '
                 '-c 0.6 -m mothur -r /baz/reference_seqs.fasta -t '
                 '/baz/id_to_taxonomy.txt')],
               [('Adding taxa (Mothur, 0.6 confidence)',
                 'add_taxa.py -i /foo/bar/otu_table.biom -o '
                 '/foo/bar/mothur_0.6/otu_table_w_taxa.biom -t '
                 '/foo/bar/mothur_0.6/rep_set_tax_assignments.txt')],
               [('Summarizing taxa (Mothur, 0.6 confidence)',
                 'summarize_taxa.py -i /foo/bar/mothur_0.6/otu_table_w_taxa.biom -o '
                 '/foo/bar/mothur_0.6')]]

        obs = _generate_mothur_commands('/foo/bar', '/foo/bar/rep_set.fna',
                '/baz/reference_seqs.fasta', '/baz/id_to_taxonomy.txt',
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
