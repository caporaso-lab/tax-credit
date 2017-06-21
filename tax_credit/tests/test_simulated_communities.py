#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2014--, tax-credit development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


from unittest import TestCase, main
from shutil import rmtree
from os import makedirs
from tempfile import mkdtemp
from os.path import join, exists
from tax_credit.simulated_communities import (generate_simulated_communities)


class EvalFrameworkTests(TestCase):

    @classmethod
    def setUpClass(cls):
        _table1 = '\n'.join([
            '#SampleID\tk__Bacteria; p__Proteobacteria; c__Gammaproteobacteria'
            '; o__Legionellales; f__Legionellaceae; g__Legionella; s__\t'
            'k__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacter'
            'iales; f__Flavobacteriaceae; g__Flavobacterium; s__',
            's1\t0.5\t0.5',
            's2\t0.1\t0.9'])

        cls.table2 = '\n'.join(['#SampleID\t0001\t0003',
                                's1\t0.5\t0.5',
                                's2\t0.1\t0.9\n'])

        _ref1 = '\n'.join([
            '0001\tk__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; '
            'o__Legionellales; f__Legionellaceae; g__Legionella; s__',
            '0003\tk__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavo'
            'bacteriales; f__Flavobacteriaceae; g__Flavobacterium; s__'])

        cls.seqs1 = '\n'.join(['>0001',
                               'ACTAGTAGTTGAC',
                               '>0003',
                               'ATCGATGCATGCA\n'])

        cls.tmpdir = mkdtemp()
        testdir = join(cls.tmpdir, 'sim_test')
        comm_dir = 'blob'
        cls.testpath = join(testdir, comm_dir)
        if not exists(cls.testpath):
            makedirs(cls.testpath)

        tab_fp = join(cls.testpath, 'expected-composition.txt')
        with open(tab_fp, 'w') as out:
            out.write(_table1)

        ref_fp = join(cls.testpath, 'ref1.tmp')
        with open(ref_fp, 'w') as out:
            out.write(_ref1)

        seqs_fp = join(cls.testpath, 'seqs1.tmp')
        with open(seqs_fp, 'w') as out:
            out.write(cls.seqs1)

        refs = {'ref1': (seqs_fp, ref_fp)}
        generate_simulated_communities(testdir, [(comm_dir, 'ref1')], refs, 1)

    def test_generate_simulated_communities(self):
        with open(join(self.testpath, 'simulated-composition.txt'), 'r') as sc:
            self.assertEqual(sc.read(), self.table2)
        with open(join(self.testpath, 'simulated-seqs.fna'), 'r') as sq:
            self.assertEqual(sq.read(), self.seqs1)

    @classmethod
    def tearDownClass(cls):
        rmtree(cls.tmpdir)


if __name__ == "__main__":
    main()
