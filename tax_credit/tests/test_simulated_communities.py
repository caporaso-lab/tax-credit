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
from os.path import join, exists
from tax_credit.simulated_communities import (generate_simulated_communities)


class EvalFrameworkTests(TestCase):

    def setUp(self):
        self.table2 = table2.split('\n')

    def test_generate_simulated_communities(self):
        testdir = 'sim_test'
        comm_dir = 'blob'
        testpath=join(testdir, comm_dir)
        if not exists(testpath):
            makedirs(testpath)

        tab_fp = join(testpath, 'expected-composition.txt')
        with open(tab_fp, 'w') as out:
            out.write(table1)

        ref_fp = join(testpath, 'ref1.tmp')
        with open(ref_fp, 'w') as out:
            out.write(ref1)

        seqs_fp = join(testpath, 'seqs1.tmp')
        with open(seqs_fp, 'w') as out:
            out.write(seqs1)

        refs = {'ref1': (seqs_fp, ref_fp)}
        generate_simulated_communities(testdir, [(comm_dir, 'ref1')], refs, 1)

        with open(join(testpath, 'simulated-composition.txt'), 'r') as comp:
            self.assertEqual(comp.read(), table2)
        with open(join(testpath, 'simulated-seqs.fna'), 'r') as sq:
            self.assertEqual(sq.read(), seqs1)

        rmtree(testdir)


table1 = """#SampleID\tk__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Legionellales; f__Legionellaceae; g__Legionella; s__\tk__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__Flavobacteriaceae; g__Flavobacterium; s__
s1\t0.5\t0.5
s2\t0.1\t0.9
"""

table2 = """#SampleID\t0001\t0003
s1\t0.5\t0.5
s2\t0.1\t0.9
"""

ref1 = """0001\tk__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Legionellales; f__Legionellaceae; g__Legionella; s__
0003\tk__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__Flavobacteriaceae; g__Flavobacterium; s__"""

seqs1 = """>0001
ACTAGTAGTTGAC
>0003
ATCGATGCATGCA
"""


if __name__ == "__main__":
    main()
