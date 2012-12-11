#!/usr/bin/env python

__author__ = "Kyle Patnode"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Kyle Patnode","Jai Ram Rideout","Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Kyle Patnode"
__email__ = "kpatnode1@gmail.com"
__status__ = "Development"

"""Test suite for the generate_taxa_compare_table.py module."""

from taxcompare.generate_taxa_compare_table import *
from os import makedirs, getcwd, chdir
from shutil import rmtree
from tempfile import mkdtemp
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files
from qiime.parse import parse_taxa_summary_table
from qiime.test import initiate_timeout, disable_timeout
from qiime.util import get_qiime_temp_dir, get_tmp_filename


class GenerateTaxaCompareTableTests(TestCase):
    """Tests for the generate_taxa_compare_table.py module."""

    def setUp(self):
        """Set up files/environment that will be used by the tests."""
        # The prefix to use for temporary files. This prefix may be added to,
        # but all temp dirs and files created by the tests will have this
        # prefix at a minimum.
        self.prefix = 'generate_taxa_compare_table_tests'

        self.start_dir = getcwd()
        self.dirs_to_remove = []
        self.files_to_remove = []
        
        self.tmp_dir = get_qiime_temp_dir()
        if not exists(self.tmp_dir):
            makedirs(self.tmp_dir)
            # if test creates the temp dir, also remove it
            self.dirs_to_remove.append(self.tmp_dir)

        # setup temporary root input directory
        self.root_dir = mkdtemp(dir=self.tmp_dir,
                                prefix='%s_root_dir_' %self.prefix)
        self.dirs_to_remove.append(self.root_dir)

        L18S_dir = '/L18S-1/blast_1.0/'
        makedirs(self.root_dir+L18S_dir)
        self.L18S_fp = self.root_dir+L18S_dir+'/otu_table_mc2_w_taxa_L5.txt'
        with open(self.L18S_fp, 'w') as f:
            f.writelines(L18S_L5_blast_one_multiple_assign_output)
        self.files_to_remove.append(self.L18S_fp)

        # setup temporary key directory
        self.key_dir = mkdtemp(dir=self.tmp_dir,
                                prefix='%s_key_dir_' %self.prefix)
        self.dirs_to_remove.append(self.key_dir)
        self.key_fp = self.key_dir+'/L18S_key.txt'
        with open(self.key_fp, 'w') as f:
            f.writelines(L18S_key)
        self.files_to_remove.append(self.key_fp)
        self.bad_key = self.key_dir+'/L18S_key.txt'

        # setup temporary output directory
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

    def test_generate_taxa_compare_table_method(self):
        """Functions correctly using standard valid input data."""
        exp = [{}, {}, {}, {'L18s-1': {'blast_1.0': ('-0.2336', '-0.7924')}}, {}]

        obs = generate_taxa_compare_table(self.root_dir, self.key_dir)

        self.assertEqual(obs, exp)

    #Test bad generate_taxa_compare_table input
    def test_invalid_generate_taxa_compare_table_input(self):
        """Test that errors are thrown using various types of invalid input."""

        out_dir = self.output_dir

        # The root directory doesn't exist.
        self.assertRaises(WorkflowError, generate_taxa_compare_table,
                '/foobarbaz', out_dir)

        # Invalid levels
        self.assertRaises(WorkflowError, generate_taxa_compare_table, out_dir, out_dir, 'foo')
        self.assertRaises(WorkflowError, generate_taxa_compare_table, out_dir, out_dir, [7])
        self.assertRaises(WorkflowError, generate_taxa_compare_table, out_dir, out_dir, [1,2,3,4,5,6])

    def test_valid_format_output(self):
        """Functions correctly using standard valid input data"""
        exp = [[], [], [], ['P,S\tblast_1.0\trdp_0.8\n', 'Broad1\tN/A\t-0.1236,-0.7477\t\n',
              'L18s-1\t-0.2336,-0.7924\tN/A\t\n'], []]

        obs = format_output([{}, {}, {}, {'L18s-1': {'blast_1.0': ('-0.2336', '-0.7924')},
              'Broad1': {'rdp_0.8': ('-0.1236', '-0.7477')}}, {}], ',')
        self.assertEqual(obs, exp)

    def test_valid_get_key_files_input(self):
        """Functions correctly using standard valid input data. Also checks to make sure get_key_files isn't grabbing backups."""
        obs = get_key_files(self.key_dir)
        # Test key and incidentally length
        self.assertTrue(obs.keys() == ['L18s'])
        # Test value
        self.assertTrue(obs.values()[0].endswith('L18S_key.txt'))

    def test_invalid_get_key_files_input(self):
        """Test that errors are thrown using various types of invalid input."""
        # Key directory does not exist
        self.assertRaises(WorkflowError, get_key_files, '/foobarbaz')

        # Key directory is empty
        self.assertRaises(WorkflowError, get_key_files, self.output_dir)

    def test_valid_get_coefficients_input(self):
        """Functions correctly using standard valid input data"""
        exp = ('-0.2336','-0.7924')

        run = parse_taxa_summary_table(open(self.L18S_fp, 'U'))
        key = parse_taxa_summary_table(open(self.key_fp, 'U'))

        obs = get_coefficients(run, key)

        self.assertEqual(obs, exp)

L18S_key = \
"""Taxon	EUK.Mock.1
Eukaryota;Fungi;Ascomycota;Saccharomycetes;Saccharomycetales;Incertae_sedis;Candida;Candida_albicans	0.083333333
Eukaryota;Fungi;Ascomycota;Saccharomycetes;Saccharomycetales;Incertae_sedis;Candida;Candida_tropicalis	0.083333333
Eukaryota;Fungi;Ascomycota;Saccharomycetes;Saccharomycetales;Saccharomycetaceae;Candida;Candida_lusitaniae	0.083333333
Eukaryota;Fungi;Ascomycota;Saccharomycetes;Saccharomycetales;Saccharomycetaceae;Saccharomyces;Saccharomyces_cerevisiae	0.083333333
Eukaryota;Fungi;Ascomycota;Sordariomycetes;Hypocreales;Nectriaceae;Fusarium;Fusarium_oxysporum	0.083333333
Eukaryota;Fungi;Ascomycota;Eurotiomycetes;Onygenales;Onygenaceae;Coccidioides;Coccidioides_immitis	0.083333333
Eukaryota;Fungi;Ascomycota;Eurotiomycetes;Onygenales;Arthrodermataceae;Trichophyton	0.083333333
Eukaryota;Fungi;Ascomycota;Schizosaccharomycetes;Schizosaccharomycetales;Schizosaccharomycetaceae;Schizosaccharomyces;Schizosaccharomyces_pombe	0.083333333
Eukaryota;Fungi;Basidiomycota;Tremellomycetes;Tremellales;Tremellaceae;Cryptococcus;Cryptococcus_neoformans	0.083333333
Eukaryota;Fungi;Incertae_sedis;Incertae_sedis;Mucorales;Incertae_sedis;Rhizopus;Rhizopus_oryzae	0.083333333
Eukaryota;Fungi;Chytridiomycota;Chytridiomycetes;Rhizophydiales;Incertae_sedis;Batrachochytrium;Batrachochytrium_dendrobatidis	0.083333333
Eukaryota;Metazoa;Chordata;Craniata;Vertebrata;Euteleostomi;Mammalia;Eutheria;Euarchontoglires;Primates;Haplorrhini;Catarrhini;Hominidae;Homo;Homo sapiens	0.083333333
"""

L18S_L5_blast_one_multiple_assign_output = \
"""Taxon	EUK.Mock.1
Bacteria;Firmicutes;Bacilli;Lactobacillales;Streptococcaceae	0.174805505122
Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae	1.88144984525e-05
Bacteria;Fusobacteria;Fusobacteria;Fusobacteriales;Fusobacteriaceae	3.7628996905e-05
Bacteria;Proteobacteria;Betaproteobacteria;Neisseriales;Neisseriaceae	0.0233958288257
Bacteria;Proteobacteria;Betaproteobacteria;SC-I-84;uncultured bacterium	1.88144984525e-05
Bacteria;Proteobacteria;Epsilonproteobacteria;Campylobacterales;Campylobacteraceae	0.00738469064261
Bacteria;Proteobacteria;Gammaproteobacteria;Aeromonadales;Aeromonadaceae	1.88144984525e-05
Bacteria;Proteobacteria;Gammaproteobacteria;Chromatiales;Chromatiaceae	1.88144984525e-05
Bacteria;Proteobacteria;Gammaproteobacteria;Legionellales;Coxiellaceae	1.88144984525e-05
Bacteria;Proteobacteria;Gammaproteobacteria;Order Incertae Sedis;Family Incertae Sedis	1.88144984525e-05
Bacteria;Proteobacteria;Gammaproteobacteria;Pasteurellales;Pasteurellaceae	0.0753896952992
Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae	0.000555027704349
Bacteria;Proteobacteria;Gammaproteobacteria;Thiotrichales;EV818SWSAP88	0.00067732194429
Bacteria;Proteobacteria;Gammaproteobacteria;Thiotrichales;Thiotrichaceae	1.88144984525e-05
Bacteria;Proteobacteria;SPOTSOCT00m83;uncultured marine bacterium;Other	1.88144984525e-05
Eukaryota;Fungi;Chytridiomycota;Chytridiomycetes;Chytridiales	0.0160581744292
Eukaryota;Fungi;Chytridiomycota;Chytridiomycetes;Cladochytriales	0.000103479741489
Eukaryota;Fungi;Dikarya;Ascomycota;Pezizomycotina	0.100808082709
Eukaryota;Fungi;Dikarya;Ascomycota;Saccharomycotina	0.305660341859
Eukaryota;Fungi;Dikarya;Ascomycota;mitosporic Ascomycota	1.88144984525e-05
Eukaryota;Fungi;Dikarya;Basidiomycota;Agaricomycotina	0.0625958363515
Eukaryota;Fungi;Dikarya;Basidiomycota;Pucciniomycotina	6.58507445838e-05
Eukaryota;Fungi;Dikarya;environmental samples;uncultured Dikarya	8.46652430363e-05
Eukaryota;Fungi;Fungi incertae sedis;Basal fungal lineages;Mucoromycotina	0.208125981882
Eukaryota;Fungi;Glomeromycota;Glomeromycetes;Glomerales	2.82217476788e-05
Eukaryota;Metazoa;Arthropoda;Myriapoda;Diplopoda	1.88144984525e-05
Eukaryota;Metazoa;Chordata;Craniata;Vertebrata	0.00895570126339
Eukaryota;Metazoa;Cnidaria;Anthozoa;Octocorallia	1.88144984525e-05
Eukaryota;Viridiplantae;Streptophyta;Embryophyta;Tracheophyta	0.000884281427268
Eukaryota;environmental samples;uncultured eukaryote;Other;Other	0.0123046819879
No blast hit;Other;Other;Other;Other	0.00187204259602
"""


if __name__ == "__main__":
    main()
