biom add-metadata#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout", "Zack Ellett", "Kyle Patnode"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"
__status__ = "Development"

"""Test suite for the multiple_assign_taxonomy.py module."""

from os import makedirs, getcwd, chdir
from os.path import exists, join
from shutil import rmtree
from tempfile import mkdtemp, NamedTemporaryFile
from cogent.util.misc import remove_files
from cogent.util.unit_test import TestCase, main
from qiime.test import initiate_timeout, disable_timeout
from qiime.util import get_qiime_temp_dir, get_tmp_filename
from qiime.workflow.util import WorkflowError

from taxcompare.multiple_assign_taxonomy import (
        assign_taxonomy_multiple_times,
        _directory_check,
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
                                  prefix='%s_output_dir_' % self.prefix)
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
        out_dir = self.output_dir

        # The output directory already exists, and we aren't in force mode.
        self.assertRaises(WorkflowError, assign_taxonomy_multiple_times,
                ['/foo', '/bar/'], out_dir, ['rdp', 'blast'],
                '/foo/ref_seqs.fasta', '/foo/tax.txt', confidences=[0.6, 0.7],
                e_values=[1e-3, 1e-10])

        # The input directories don't exist.
        self.assertRaises(WorkflowError, assign_taxonomy_multiple_times,
                ['/foobarbaz', '/foobarbaz2/'], out_dir, ['rdp', 'blast'],
                '/foo/ref_seqs.fasta', '/foo/tax.txt', confidences=[0.6, 0.7],
                e_values=[1e-3, 1e-10], force=True)

        # Invalid assignment method.
        self.assertRaises(WorkflowError, assign_taxonomy_multiple_times,
                          [out_dir, out_dir], out_dir, ['foo', 'rdp'],
                          '/foo/ref_seqs.fasta', '/foo/tax.txt',
                          confidences=[0.6, 0.7], e_values=[1e-3, 1e-10],
                          force=True)

        # Not enough custom regexes.
        self.assertRaises(WorkflowError, assign_taxonomy_multiple_times,
                          [out_dir, out_dir], out_dir, ['rdp'],
                          '/foo/ref_seqs.fasta', '/foo/tax.txt',
                          confidences=[0.6, 0.7], e_values=[1e-3, 1e-10],
                          force=True, rtax_amplicon_id_regexes=['foo*'])

    def test_assign_taxonomy_multiple_times_invalid_input_rdp(self):
        """Test that errors are thrown using invalid input for RDP."""
        out_dir = self.output_dir

        # RDP confidences are missing.
        self.assertRaises(WorkflowError, assign_taxonomy_multiple_times,
                          [out_dir, out_dir], out_dir, ['rdp', 'rdp'],
                          '/foo/ref_seqs.fasta', '/foo/tax.txt',
                          e_values=[1e-3, 1], force=True)

    def test_assign_taxonomy_multiple_times_invalid_input_blast(self):
        """Test that errors are thrown using invalid input for BLAST."""
        out_dir = self.output_dir

        # BLAST E-values are missing.
        self.assertRaises(WorkflowError, assign_taxonomy_multiple_times,
                          [out_dir, out_dir], out_dir, ['blast', 'blast'],
                          '/foo/ref_seqs.fasta', '/foo/tax.txt',
                          confidences=[0.88], force=True)

    def test_assign_taxonomy_multiple_times_invalid_input_mothur(self):
        """Test that errors are thrown using invalid input for Mothur."""
        out_dir = self.output_dir

        # Mothur confidences are missing.
        self.assertRaises(WorkflowError, assign_taxonomy_multiple_times,
                          [out_dir, out_dir], out_dir, ['mothur', 'mothur'],
                          '/foo/ref_seqs.fasta', '/foo/tax.txt',
                          rtax_modes=['single', 'paired'], force=True)

    def test_assign_taxonomy_multiple_times_invalid_input_rtax(self):
        """Test that errors are thrown using invalid input for RTAX."""
        out_dir = self.output_dir

        # RTAX modes are missing.
        self.assertRaises(WorkflowError, assign_taxonomy_multiple_times,
                          [out_dir, out_dir], out_dir, ['rtax', 'rtax'],
                          '/foo/ref_seqs/fasta', '/foo/tax.txt',
                          force=True)

        # RTAX modes are invalid.
        self.assertRaises(WorkflowError, assign_taxonomy_multiple_times,
                          [out_dir, out_dir], out_dir, ['rtax', 'rtax'],
                          '/foo/ref_seqs/fasta', '/foo/tax.txt',
                          rtax_modes=['single', 'piared'], force=True)

    def test_directory_check(self):
        """Test that directory names are generated properly."""
        exp = ("output_dir/method_X","output_dir/method_X.tmp")

        obs = _directory_check("output_dir", "method_", "X")
        self.assertEqual(obs, exp)

    def test_generate_rdp_commands(self):
        """Functions correctly using standard valid input data."""
        exp = [[('Assigning taxonomy (RDP, confidence: 0.8)',
                 'assign_taxonomy.py -i /foo/bar/rep_set.fna -o /foo/bar/rdp_0.8.tmp '
                 '-c 0.8 -m rdp -r /baz/reference_seqs.fasta -t /baz/id_to_taxonomy.txt')],
                 [('Adding metadata (RDP, confidence: 0.8)',
                 'biom add-metadata -i /foo/bar/otu_table.biom -o '
                 '/foo/bar/rdp_0.8.tmp/otu_table_w_taxa.biom '
                 '--observation-metadata-fp '
                 '/foo/bar/rdp_0.8.tmp/rep_set_tax_assignments.txt '
                 '--sc-separated taxonomy --observation-header OTUID,taxonomy')],
                 [('Summarizing taxa (RDP, confidence: 0.8)',
                 'summarize_taxa.py -i /foo/bar/rdp_0.8.tmp/otu_table_w_taxa.biom -o '
                 '/foo/bar/rdp_0.8.tmp')],
                 [('Renaming output directory (RDP, confidence: 0.8)',
                 'mv /foo/bar/rdp_0.8.tmp /foo/bar/rdp_0.8')],
                 [('Assigning taxonomy (RDP, confidence: 0.6)',
                 'assign_taxonomy.py -i /foo/bar/rep_set.fna -o /foo/bar/rdp_0.6.tmp '
                 '-c 0.6 -m rdp -r /baz/reference_seqs.fasta -t /baz/id_to_taxonomy.txt')],
                 [('Adding metadata (RDP, confidence: 0.6)',
                 'biom add-metadata -i /foo/bar/otu_table.biom -o '
                 '/foo/bar/rdp_0.6.tmp/otu_table_w_taxa.biom '
                 '--observation-metadata-fp '
                 '/foo/bar/rdp_0.6.tmp/rep_set_tax_assignments.txt '
                 '--sc-separated taxonomy --observation-header OTUID,taxonomy')],
                 [('Summarizing taxa (RDP, confidence: 0.6)',
                 'summarize_taxa.py -i /foo/bar/rdp_0.6.tmp/otu_table_w_taxa.biom -o '
                 '/foo/bar/rdp_0.6.tmp')],
                 [('Renaming output directory (RDP, confidence: 0.6)',
                 'mv /foo/bar/rdp_0.6.tmp /foo/bar/rdp_0.6')]]

        obs = _generate_rdp_commands('/foo/bar', '/foo/bar/rep_set.fna',
                '/baz/reference_seqs.fasta', '/baz/id_to_taxonomy.txt',
                '/foo/bar/otu_table.biom', [0.80, 0.60])
        self.assertEqual(obs, exp)

        # Test rdp_max_memory.
        exp = [[('Assigning taxonomy (RDP, confidence: 0.8)',
                 'assign_taxonomy.py -i /foo/bar/rep_set.fna -o /foo/bar/rdp_0.8.tmp '
                 '-c 0.8 -m rdp -r /baz/reference_seqs.fasta -t '
                 '/baz/id_to_taxonomy.txt --rdp_max_memory 2')],
                 [('Adding metadata (RDP, confidence: 0.8)',
                 'biom add-metadata -i /foo/bar/otu_table.biom -o '
                 '/foo/bar/rdp_0.8.tmp/otu_table_w_taxa.biom '
                 '--observation-metadata-fp '
                 '/foo/bar/rdp_0.8.tmp/rep_set_tax_assignments.txt '
                 '--sc-separated taxonomy --observation-header OTUID,taxonomy')],
                 [('Summarizing taxa (RDP, confidence: 0.8)',
                 'summarize_taxa.py -i /foo/bar/rdp_0.8.tmp/otu_table_w_taxa.biom -o '
                 '/foo/bar/rdp_0.8.tmp')],
                 [('Renaming output directory (RDP, confidence: 0.8)',
                 'mv /foo/bar/rdp_0.8.tmp /foo/bar/rdp_0.8')]]

        obs = _generate_rdp_commands('/foo/bar', '/foo/bar/rep_set.fna',
                '/baz/reference_seqs.fasta', '/baz/id_to_taxonomy.txt',
                '/foo/bar/otu_table.biom', [0.80], rdp_max_memory=2)
        self.assertEqual(obs, exp)

        # Test skips directory that already exists.
        try:
            makedirs(join(self.output_dir, 'rdp_0.8'))
        except OSError:
            pass

        obs = _generate_rdp_commands(self.output_dir, '/foo/bar/rep_set.fna',
                '/baz/reference_seqs.fasta', '/baz/id_to_taxonomy.txt',
                '/foo/bar/otu_table.biom', [0.80])
        self.assertEqual(obs, [])

    def test_generate_blast_commands(self):
        """Functions correctly using standard valid input data."""
        exp = [[('Assigning taxonomy (BLAST, E: 0.002)',
                 'assign_taxonomy.py -i /foo/bar/rep_set.fna -o /foo/bar/blast_0.002.tmp -e 0.002 '
                 '-m blast -r /baz/reference_seqs.fasta -t /baz/id_to_taxonomy.txt')],
                 [('Adding metadata (BLAST, E: 0.002)',
                 'biom add-metadata -i /foo/bar/otu_table.biom -o '
                 '/foo/bar/blast_0.002.tmp/otu_table_w_taxa.biom '
                 '--observation-metadata-fp '
                 '/foo/bar/blast_0.002.tmp/rep_set_tax_assignments.txt '
                 '--sc-separated taxonomy --observation-header OTUID,taxonomy')],
                 [('Summarizing taxa (BLAST, E: 0.002)',
                 'summarize_taxa.py -i /foo/bar/blast_0.002.tmp/otu_table_w_taxa.biom '
                 '-o /foo/bar/blast_0.002.tmp')],
                 [('Renaming output directory (BLAST, E: 0.002)',
                 'mv /foo/bar/blast_0.002.tmp /foo/bar/blast_0.002')],
                 [('Assigning taxonomy (BLAST, E: 0.005)',
                 'assign_taxonomy.py -i /foo/bar/rep_set.fna -o /foo/bar/blast_0.005.tmp -e 0.005 '
                 '-m blast -r /baz/reference_seqs.fasta -t /baz/id_to_taxonomy.txt')],
                 [('Adding metadata (BLAST, E: 0.005)',
                 'biom add-metadata -i /foo/bar/otu_table.biom -o '
                 '/foo/bar/blast_0.005.tmp/otu_table_w_taxa.biom '
                 '--observation-metadata-fp '
                 '/foo/bar/blast_0.005.tmp/rep_set_tax_assignments.txt '
                 '--sc-separated taxonomy --observation-header OTUID,taxonomy')],
                 [('Summarizing taxa (BLAST, E: 0.005)',
                 'summarize_taxa.py -i /foo/bar/blast_0.005.tmp/otu_table_w_taxa.biom '
                 '-o /foo/bar/blast_0.005.tmp')],
                 [('Renaming output directory (BLAST, E: 0.005)',
                 'mv /foo/bar/blast_0.005.tmp /foo/bar/blast_0.005')]]

        obs = _generate_blast_commands('/foo/bar', '/foo/bar/rep_set.fna',
                '/baz/reference_seqs.fasta', '/baz/id_to_taxonomy.txt',
                '/foo/bar/otu_table.biom', [0.002, 0.005])
        self.assertEqual(obs, exp)

    def test_generate_mothur_commands(self):
        """Functions correctly using standard valid input data."""
        exp = [[('Assigning taxonomy (Mothur, confidence: 0.8)',
                 'assign_taxonomy.py -i /foo/bar/rep_set.fna -o /foo/bar/mothur_0.8.tmp '
                 '-c 0.8 -m mothur -r /baz/reference_seqs.fasta -t '
                 '/baz/id_to_taxonomy.txt')],
                 [('Adding metadata (Mothur, confidence: 0.8)',
                 'biom add-metadata -i /foo/bar/otu_table.biom -o '
                 '/foo/bar/mothur_0.8.tmp/otu_table_w_taxa.biom '
                 '--observation-metadata-fp '
                 '/foo/bar/mothur_0.8.tmp/rep_set_tax_assignments.txt '
                 '--sc-separated taxonomy --observation-header OTUID,taxonomy')],
                 [('Summarizing taxa (Mothur, confidence: 0.8)',
                 'summarize_taxa.py -i /foo/bar/mothur_0.8.tmp/otu_table_w_taxa.biom -o '
                 '/foo/bar/mothur_0.8.tmp')],
                 [('Renaming output directory (Mothur, confidence: 0.8)',
                 'mv /foo/bar/mothur_0.8.tmp /foo/bar/mothur_0.8')],
                 [('Assigning taxonomy (Mothur, confidence: 0.6)',
                 'assign_taxonomy.py -i /foo/bar/rep_set.fna -o /foo/bar/mothur_0.6.tmp '
                 '-c 0.6 -m mothur -r /baz/reference_seqs.fasta -t '
                 '/baz/id_to_taxonomy.txt')],
                 [('Adding metadata (Mothur, confidence: 0.6)',
                 'biom add-metadata -i /foo/bar/otu_table.biom -o '
                 '/foo/bar/mothur_0.6.tmp/otu_table_w_taxa.biom '
                 '--observation-metadata-fp '
                 '/foo/bar/mothur_0.6.tmp/rep_set_tax_assignments.txt '
                 '--sc-separated taxonomy --observation-header OTUID,taxonomy')],
                 [('Summarizing taxa (Mothur, confidence: 0.6)',
                 'summarize_taxa.py -i /foo/bar/mothur_0.6.tmp/otu_table_w_taxa.biom -o '
                 '/foo/bar/mothur_0.6.tmp')],
                 [('Renaming output directory (Mothur, confidence: 0.6)',
                 'mv /foo/bar/mothur_0.6.tmp /foo/bar/mothur_0.6')]]

        obs = _generate_mothur_commands('/foo/bar', '/foo/bar/rep_set.fna',
                '/baz/reference_seqs.fasta', '/baz/id_to_taxonomy.txt',
                '/foo/bar/otu_table.biom', [0.80, 0.60])
        self.assertEqual(obs, exp)

    def test_generate_rtax_commands(self):
        """Functions correctly using standard valid input data."""
        exp = [[('Assigning taxonomy (RTAX, mode: single)',
                 'assign_taxonomy.py -i /foo/bar/rep_set.fna -o /foo/bar/rtax_single.tmp '
                 '-m rtax -r /baz/reference_seqs.fasta -t /baz/id_to_taxonomy.txt '
                 '--read_1_seqs_fp /foo/bar/read_1_seqs.fna')],
                 [('Adding metadata (RTAX, mode: single)',
                 'biom add-metadata -i /foo/bar/otu_table.biom -o '
                 '/foo/bar/rtax_single.tmp/otu_table_w_taxa.biom '
                 '--observation-metadata-fp '
                 '/foo/bar/rtax_single.tmp/rep_set_tax_assignments.txt '
                 '--sc-separated taxonomy --observation-header OTUID,taxonomy')],
                 [('Summarizing taxa (RTAX, mode: single)',
                 'summarize_taxa.py -i /foo/bar/rtax_single.tmp/otu_table_w_taxa.biom '
                 '-o /foo/bar/rtax_single.tmp')],
                 [('Renaming output directory (RTAX, mode: single)',
                 'mv /foo/bar/rtax_single.tmp /foo/bar/rtax_single')],
                 [('Assigning taxonomy (RTAX, mode: paired)',
                 'assign_taxonomy.py -i /foo/bar/rep_set.fna -o /foo/bar/rtax_paired.tmp '
                 '-m rtax -r /baz/reference_seqs.fasta -t /baz/id_to_taxonomy.txt '
                 '--read_1_seqs_fp /foo/bar/read_1_seqs.fna '
                 '--read_2_seqs_fp /foo/bar/read_2_seqs.fna --single_ok')],
                 [('Adding metadata (RTAX, mode: paired)',
                 'biom add-metadata -i /foo/bar/otu_table.biom -o '
                 '/foo/bar/rtax_paired.tmp/otu_table_w_taxa.biom '
                 '--observation-metadata-fp '
                 '/foo/bar/rtax_paired.tmp/rep_set_tax_assignments.txt '
                 '--sc-separated taxonomy --observation-header OTUID,taxonomy')],
                 [('Summarizing taxa (RTAX, mode: paired)',
                 'summarize_taxa.py -i /foo/bar/rtax_paired.tmp/otu_table_w_taxa.biom '
                 '-o /foo/bar/rtax_paired.tmp')],
                 [('Renaming output directory (RTAX, mode: paired)',
                 'mv /foo/bar/rtax_paired.tmp /foo/bar/rtax_paired')]]

        obs = _generate_rtax_commands('/foo/bar', '/foo/bar/rep_set.fna',
                '/baz/reference_seqs.fasta', '/baz/id_to_taxonomy.txt',
                '/foo/bar/otu_table.biom', ['single', 'paired'],
                '/foo/bar/read_1_seqs.fna', '/foo/bar/read_2_seqs.fna')
        self.assertEqual(obs, exp)

    def test_generate_rtax_commands_custom_regexes(self):
        """Correctly passes custom regexes through to RTAX."""
        exp = [[('Assigning taxonomy (RTAX, mode: single)',
                 'assign_taxonomy.py -i /foo/bar/rep_set.fna -o /foo/bar/rtax_single.tmp -m rtax -r /baz/reference_seqs.fasta -t /baz/id_to_taxonomy.txt --read_1_seqs_fp /foo/bar/read_1_seqs.fna --read_id_regex \'f*\' --amplicon_id_regex \'*b\' --header_id_regex \'b(*)z\'')],
               [('Assigning taxonomy (RTAX, mode: paired)',
                 'assign_taxonomy.py -i /foo/bar/rep_set.fna -o /foo/bar/rtax_paired.tmp ' '-m rtax -r /baz/reference_seqs.fasta -t /baz/id_to_taxonomy.txt ' '--read_1_seqs_fp /foo/bar/read_1_seqs.fna ' '--read_2_seqs_fp /foo/bar/read_2_seqs.fna --read_id_regex \'f*\' --amplicon_id_regex \'*b\' --header_id_regex \'b(*)z\' --single_ok')]]

        obs = _generate_rtax_commands('/foo/bar', '/foo/bar/rep_set.fna',
                '/baz/reference_seqs.fasta', '/baz/id_to_taxonomy.txt',
                '/foo/bar/otu_table.biom', ['single', 'paired'],
                '/foo/bar/read_1_seqs.fna', '/foo/bar/read_2_seqs.fna', 'f*',
                '*b', 'b(*)z')

        self.assertEqual(len(obs), 8)
        self.assertEqual(obs[0], exp[0])
        self.assertEqual(obs[4], exp[1])

    def test_generate_taxa_processing_commands(self):
        """Functions correctly using standard valid input data."""
        exp = ([('Adding metadata (RDP, confidence: 0.8)',
                 'biom add-metadata -i /foo/otu_table.biom -o '
                 '/foo/rdp_0.8/otu_table_w_taxa.biom '
                 '--observation-metadata-fp '
                 '/foo/rdp_0.8/rep_set_tax_assignments.txt '
                 '--sc-separated taxonomy --observation-header OTUID,taxonomy')],
               [('Summarizing taxa (RDP, confidence: 0.8)',
                 'summarize_taxa.py -i /foo/rdp_0.8/otu_table_w_taxa.biom -o '
                 '/foo/rdp_0.8')])

        obs = _generate_taxa_processing_commands('/foo/rdp_0.8',
                '/foo/rep_set.fna', '/foo/otu_table.biom',
                'RDP, confidence: 0.8')
        self.assertEqual(obs, exp)


if __name__ == "__main__":
    main()
