#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2014--, tax-credit development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


from unittest import TestCase, main
from tempfile import mkdtemp
from os.path import join
from shutil import rmtree
from tax_credit.taxa_manipulator import (string_search,
                                         unique_lines,
                                         trim_taxonomy_strings,
                                         branching_taxa,
                                         extract_taxa_names,
                                         compile_reference_taxa,
                                         extract_rownames,
                                         extract_fasta_ids,
                                         filter_sequences,
                                         stratify_taxonomy_subsets,
                                         accept_list_or_file)


class EvalFrameworkTests(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.op11 = \
            '203525\tk__Bacteria; p__OP11; c__OP11-1; o__; f__; g__; s__'
        cls.op11_t = 'k__Bacteria; p__OP11; c__OP11-1; o__; f__; g__; s__'

        cls.table1 = [
            '229854\tk__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; '
            'o__Legionellales; f__Legionellaceae; g__Legionella; s__',
            '367523\tk__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flav'
            'obacteriales; f__Flavobacteriaceae; g__Flavobacterium; s__',
            '239330\tk__Bacteria; p__Proteobacteria; c__Deltaproteobacteria; '
            'o__Desulfuromonadales; f__Geobacteraceae; g__Geobacter; s__',
            '203525\tk__Bacteria; p__OP11; c__OP11-1; o__; f__; g__; s__']

        cls.table2 = [
            '229854\tk__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; '
            'o__Legionellales; f__Legionellaceae; g__Legionella; s__',
            '229854\tk__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; '
            'o__Legionellales; f__Legionellaceae; g__Legionella; s__',
            '367523\tk__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flav'
            'obacteriales; f__Flavobacteriaceae; g__Flavobacterium; s__',
            '367523\tk__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flav'
            'obacteriales; f__Flavobacteriaceae; g__Flavobacterium; s__',
            '239330\tk__Bacteria; p__Proteobacteria; c__Deltaproteobacteria; '
            'o__Desulfuromonadales; f__Geobacteraceae; g__Geobacter; s__',
            '203525\tk__Bacteria; p__OP11; c__OP11-1; o__; f__; g__; s__',
            '239330\tk__Bacteria; p__Proteobacteria; c__Deltaproteobacteria; '
            'o__Desulfuromonadales; f__Geobacteraceae; g__Geobacter; s__']

        cls.table3 = [
            '229854\tk__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; '
            'o__Legionellales; f__Legionellaceae; g__Legionella; s__',
            '367523\tk__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flav'
            'obacteriales; f__Flavobacteriaceae; g__Flavobacterium; s__',
            '239330\tk__Bacteria; p__Proteobacteria; c__Deltaproteobacteria; '
            'o__Desulfuromonadales; f__Geobacteraceae; g__Geobacter; s__']

        cls.seqs1 = '\n'.join(['>229854',
                               'ACTAGTAGTTGAC',
                               '>367523',
                               'ATCGATGCATGCA',
                               '>239330',
                               'TGTGTGCTGGTAGTTAC',
                               '>203525',
                               'TGTATGCTGATGC\n'])

        cls.seqs2 = '\n'.join(['>229854',
                               'ACTAGTAGTTGAC',
                               '>367523',
                               'ATCGATGCATGCA',
                               '>239330',
                               'TGTGTGCTGGTAGTTAC\n'])

        cls.tmpdir = mkdtemp()

    def test_string_search(self):
        self.assertEqual(string_search(self.table1, '203525'), [self.op11])
        self.assertEqual(string_search(self.table1, '229854|367523|239330',
                                       discard=True), [self.op11])
        self.assertEqual(string_search(self.table1, 'OP11-1'), [self.op11])
        self.assertEqual(string_search(self.table1, 'OP11-1'), [self.op11])

    def test_unique_lines(self):
        self.assertEqual(unique_lines(self.table2), self.table1)
        self.assertEqual(unique_lines(self.table2, 'u'), [self.op11])
        self.assertEqual(unique_lines(self.table2, 'u', field=1,
                                      printfield=True), [self.op11_t])
        self.assertEqual(set(unique_lines(self.table2, 'd')), set(self.table3))

    def test_trim_taxonomy_strings(self):
        self.assertEqual(trim_taxonomy_strings(
            [self.op11], level=6), [self.op11])
        self.assertEqual(trim_taxonomy_strings(
            [self.op11], level=0), ['203525\tk__Bacteria'])

    def test_branching_taxa(self):
        self.assertEqual(branching_taxa(self.table1, field=6), [])
        self.assertEqual(set(branching_taxa(self.table1, field=1)),
                         set(self.table1))
        self.assertEqual(set(branching_taxa(self.table1, field=2)),
                         set(['229854\tk__Bacteria; p__Proteobacteria; '
                              'c__Gammaproteobacteria; o__Legionellales; '
                              'f__Legionellaceae; g__Legionella; s__',
                              '239330\tk__Bacteria; p__Proteobacteria; '
                              'c__Deltaproteobacteria; o__Desulfuromonadales; '
                              'f__Geobacteraceae; g__Geobacter; s__']))

    def test_extract_taxa_names(self):
        self.assertEqual(extract_taxa_names(
            self.table1), ['s__', 's__', 's__', 's__'])
        self.assertEqual(extract_taxa_names(
            self.table1, level=slice(0, 1), field=1),
            ['k__Bacteria', 'k__Bacteria', 'k__Bacteria', 'k__Bacteria'])

    def test_compile_reference_taxa(self):
        self.assertEqual(compile_reference_taxa(
            [self.op11]), {self.op11_t: ['203525']})
        self.assertEqual(compile_reference_taxa(
            self.table1)[self.op11_t], ['203525'])

    def test_extract_rownames(self):
        self.assertEqual(extract_rownames(self.table1),
                         {'229854', '367523', '239330', '203525'})

    def test_extract_fasta_ids(self):
        with open(join(self.tmpdir, 'seqs1.tmp'), 'w') as out:
            out.write(self.seqs1)
        self.assertEqual(extract_fasta_ids(join(self.tmpdir, 'seqs1.tmp')),
                         {'229854', '367523', '239330', '203525'})

    def test_filter_sequences(self):
        filter_sequences(join(self.tmpdir, 'seqs1.tmp'),
                         join(self.tmpdir, 'seqs.tmp'), self.table3)
        with open(join(self.tmpdir, 'seqs.tmp'), 'r') as sq:
            self.assertEqual(sq.read(), self.seqs2)
        filter_sequences(join(self.tmpdir, 'seqs1.tmp'),
                         join(self.tmpdir, 'seqs.tmp'),
                         [self.op11], keep=False)
        with open(join(self.tmpdir, 'seqs.tmp'), 'r') as sq:
            self.assertEqual(sq.read(), self.seqs2)

    def test_stratify_taxonomy_subsets(self):
        stratify_taxonomy_subsets(self.table2, 2, self.tmpdir, 'test', level=5)
        iter0 = accept_list_or_file(
            join(self.tmpdir, 'test-iter0/query_taxa.tsv'))
        iter1 = accept_list_or_file(
            join(self.tmpdir, 'test-iter1/query_taxa.tsv'))
        self.assertEqual(set(iter0) & set(iter1), set(self.table3))

    @classmethod
    def tearDownClass(cls):
        rmtree(cls.tmpdir)


if __name__ == "__main__":
    main()
