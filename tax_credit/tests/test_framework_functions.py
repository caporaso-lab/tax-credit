#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2014--, tax-credit development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os import makedirs
from os.path import exists, join
from collections import Counter

from unittest import TestCase, main
import pandas as pd
from shutil import rmtree
from tempfile import mkdtemp
from tax_credit.framework_functions import (
    generate_simulated_datasets,
    find_last_common_ancestor,
    novel_taxa_classification_evaluation,
    extract_per_level_accuracy)
from tax_credit.taxa_manipulator import (import_to_list,
                                         import_taxonomy_to_dict,
                                         extract_taxa_names)


class EvalFrameworkTests(TestCase):

    def test_generate_simulated_datasets(self):
        generate_simulated_datasets(
            self.ref_data, self.tmpdir, 100, 2, range(6, 5, -1))
        # cross-validated q0 should be ref1 and vice versa
        q0 = import_to_list(
            join(self.cvdir, 'B1-REF-iter0/query_taxa.tsv'), field=1)
        q1 = import_to_list(
            join(self.cvdir, 'B1-REF-iter1/query_taxa.tsv'), field=1)
        ref0 = import_to_list(
            join(self.cvdir, 'B1-REF-iter0/ref_taxa.tsv'), field=1)
        ref1 = import_to_list(
            join(self.cvdir, 'B1-REF-iter1/ref_taxa.tsv'), field=1)
        for q_taxon in q0:
            for r_taxon in ref1:
                if r_taxon.startswith(q_taxon):
                    break
            else:
                self.fail(q_taxon + ' not in reference taxa')
        for q_taxon in q1:
            for r_taxon in ref0:
                if r_taxon.startswith(q_taxon):
                    break
            else:
                self.fail(q_taxon + ' not in reference taxa')
        # ref0 + ref1 should equal ref
        ref = import_to_list(join(self.tmpdir, self.query_fp), field=1)
        self.assertEqual(Counter(ref), Counter(ref0 + ref1))
        # confirm that query seq IDs are not in ref for cross-val
        for i in [0, 1]:
            # test that cross-validated queries have pair in ref,
            # but keys do not match
            query_ids = import_to_list(
                join(self.cvdir, 'B1-REF-iter{0}/query_taxa.tsv'.format(i)),
                field=0)
            ref_ids = import_to_list(
                join(self.cvdir, 'B1-REF-iter{0}/ref_taxa.tsv'.format(i)),
                field=0)
            # seq ID keys do not match
            self.assertEqual(set(query_ids).intersection(set(ref_ids)), set())
            # test that novel-taxa queries have no match in ref
            query_taxa = import_taxonomy_to_dict(
                join(self.ntdir, 'B1-REF-L6-iter{0}/query_taxa.tsv'.format(i)))
            ref_taxa = import_taxonomy_to_dict(
                join(self.ntdir, 'B1-REF-L6-iter{0}/ref_taxa.tsv'.format(i)))
            for key, value in query_taxa.items():
                self.assertNotIn(key, ref_taxa)
                self.assertNotIn(value, ref_taxa.values())

    def test_find_last_common_ancestor(self):
        taxa1 = import_to_list(self.query_fp)
        taxa2 = import_to_list(self.obs_taxa_fp)
        t1 = extract_taxa_names(taxa1, level=slice(None), field=1)
        t2 = extract_taxa_names(taxa2, level=slice(None), field=1)
        for i, n in zip(range(8), [7, 6, 5, 4, 7, 3, 6, 6]):
            self.assertEqual(find_last_common_ancestor(t2[i], t1[i]), n)

    def test_novel_taxa_classification_evaluation(self):
        # test novel taxa evaluation
        results = novel_taxa_classification_evaluation(
            [self.paramdir], self.tmpdir, join(self.tmpdir, 'summary.txt'),
            test_type='novel-taxa')
        self.assertEqual(results.iloc[0]['Dataset'], 'B1-REF')
        self.assertEqual(int(results.iloc[0]['level']), 6)
        self.assertEqual(int(results.iloc[0]['iteration']), 0)
        self.assertEqual(results.iloc[0]['Method'], 'method1')
        self.assertEqual(results.iloc[0]['Parameters'], 'param1')
        self.assertAlmostEqual(results.iloc[0]['match_ratio'], 0.25)
        self.assertAlmostEqual(results.iloc[0]['overclassification_ratio'], 0.)
        self.assertAlmostEqual(
            results.iloc[0]['underclassification_ratio'], 0.625)
        self.assertAlmostEqual(
            results.iloc[0]['misclassification_ratio'], 0.125)
        self.assertAlmostEqual(results.iloc[0]['Precision'], 0.666666667)
        self.assertAlmostEqual(results.iloc[0]['Recall'], 0.25)
        self.assertAlmostEqual(
            results.iloc[0]['F-measure'], 0.363636364)
        self.assertEqual(results.iloc[0]['mismatch_level_list'], self.exp_m)
        # test cross-validated evaluation
        results = novel_taxa_classification_evaluation(
            [self.paramdir], self.tmpdir, join(self.tmpdir, 'summary.txt'),
            test_type='cross-validated')
        self.assertEqual(results.iloc[0]['Dataset'], 'B1-REF-L6')
        self.assertEqual(int(results.iloc[0]['iteration']), 0)
        self.assertEqual(results.iloc[0]['Method'], 'method1')
        self.assertEqual(results.iloc[0]['Parameters'], 'param1')
        self.assertAlmostEqual(results.iloc[0]['match_ratio'], 0.25)
        self.assertEqual(results.iloc[0]['overclassification_ratio'], 0)
        self.assertAlmostEqual(
            results.iloc[0]['underclassification_ratio'], 0.625)
        self.assertAlmostEqual(
            results.iloc[0]['misclassification_ratio'], 0.125)
        self.assertEqual(results.iloc[0]['Precision'], self.exp_p)
        self.assertEqual(results.iloc[0]['Recall'], self.exp_r)
        self.assertEqual(results.iloc[0]['F-measure'], self.exp_f)
        self.assertEqual(results.iloc[0]['mismatch_level_list'], self.exp_m)

    def test_extract_per_level_accuracy(self):
        results = novel_taxa_classification_evaluation(
            [self.paramdir], self.tmpdir, join(self.tmpdir, 'summary.txt'),
            test_type='cross-validated')
        pla = extract_per_level_accuracy(results)
        # confirm that method/dataset data propagate properly
        self.assertEqual(pla['Dataset'].unique(), ['B1-REF-L6'])
        self.assertEqual(pla['iteration'].unique(), ['0'])
        self.assertEqual(pla['Method'].unique(), ['method1'])
        self.assertEqual(pla['Parameters'].unique(), ['param1'])
        # confirm that level and per-level results are correct
        for i in range(6):
            self.assertEqual(pla.iloc[i]['level'], i + 1)
            self.assertEqual(pla.iloc[i]['Precision'], self.exp_p[i + 1])
            self.assertEqual(pla.iloc[i]['Recall'], self.exp_r[i + 1])
            self.assertEqual(pla.iloc[i]['F-measure'], self.exp_f[i + 1])
            self.assertEqual(pla.iloc[i]['match_ratio'], self.exp_r[i + 1])

    @classmethod
    def setUpClass(self):
        _ref1 = '\n'.join([
            '>179419',
            'TGAGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCATGCCTAATACATGCAAGTCGAACGAG'
            'CTTCCGTTGAATGACGTGCTTGCACTGATTTCAACAATGAAGCGAGTGGCGAACTGGTGAGTAAC'
            'ACGTGGGGAATCTGCCCAGAAGCAGGGGATAACACTTGGAAACAGGTGCTAATACCGTATAACAA'
            'CAAAATCCGCATGGATCTTGTTTGAAAGGTGGCTTCGGCTATCACTTCTGGATGATCCCGCGGCG'
            'TATTAGTTAGTTGGTGAGGTAAAGGCCCACCAAGACGATGATACGTAGCCGACCTGAGAGGGTAA'
            'TCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCA'
            'CAATGGACGAAAGTCTGATGGAGCAATGCCGCGTGAGTGAAGAAGGGTTTCGGCTCGTAAAACTC'
            'TGTTGTTAAAGAAGAACACCTTTGAGAGTAACTGTTCAAGGGTTGACGGTATTTAACCAGAAAGC'
            'CACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTG'
            'GGCGTAAAGCGAGCGCAGGCGGTTTTTTAAGTCTGATGTGAAAGCCTTCGGCTTAACCGGAGAAG'
            'TGCATCGGAAACTGGGAGACTTGAGTGCAGAAGAGGACAGTGGAACTCCATGTGTAGCGGTGGAA'
            'TGCGTAGATATATGGAAGAACACCAGTGGCGAAGGCGGCTGTCTAGTCTGTAACTGACGCTGAGG'
            'CTCGAAAGCATGGGTAGCGAACAGGATTAGATACCCTGGTAGTCCATGCCGTAAACGATGAGTGC'
            'TAAGTGTTGGAGGGTTTCCGCCCTTCAGTGCTGCAGCTAACGCATTAAGCACTCCGCCTGGGGAG'
            'TACGACCGCAAGGTTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGT'
            'TTAATTCGAAGCTACGCGAAGAACCTTACCAGGTCTTGACATCTTCTGCCAATCTTAGAGATAAG'
            'ACGTTCCCTTCGGGGACAGAATGACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTCGTGAGATGT'
            'TGGGTTAAGTCCCGCAACGAGCGCAACCCTTATTATCAGTTGCCAGCATTCAGTTGGGCACTCTG'
            'GTGAGACTGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAATCATCATGCCCCTTATGA'
            'CCTGGGCTACACACGTGCTACAATGGACGGTACAACGAGTCGCGAAGTCGTGAGGCTAAGCTAAT'
            'CTCTTAAAGCCGTTCTCAGTTCGGATTGTAGGCTGCAACTCGCCTACATGAAGTTGGAATCGCTA'
            'GTAATCGCGGATCAGCATGCCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACAC'
            'CATGAGAGTTTGTAACACCCAAAGCCGGTGAGATAACCTTCGGGAGTCAGCCGTCTAAGGTGGGA'
            'CAGATGATTAGGGTGAAGTCGTAACAAGGTAGCCGTAGGAGAACCTGCGGCTGGATCACCTCCTT'
            'TCT',
            '>1117026',
            'GAGTGGCGAACTGGTGAGTAACACGTGGGAAATCTGCCCAGAAGCAGGGGATAACACTTGGAAAC'
            'AGGTGCTAATACCGTATAACAACAAGAACCGCATGGTTCTTGTTTGAAAGGTGGTTTCGGCTATC'
            'ACTTCTGGATGATCCCGCGGCGTATTAGTTAGTTGGTGAGGTAAAGGCCCACCAAGACAATGATA'
            'CGTAGCCGACCTGAGAGGGTAATCGGCCACATTGGGACTGAGACACGGCCCAAACTCCTACGGGA'
            'GGCAGCAGTAGGGAATCTTCCACAATGGACGAAAGTCTGATGGAGCAATGCCGCGTGAGTGAAGA'
            'AGGGTTTCGGCTCGTAAAACTCTGTTGTTAAAGAAGAACACCTCTGAGAGTAACTGTTCAGGGGT'
            'TGACGGTATTTAACCAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGG'
            'CAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTTCTTAAGTCTGATGTGAAA'
            'GCCTTCGGCTTAACCGGAGAAGTGCATCGGAAACTGGGTAACTTGAGTGCAGAAGAGGACAGTGG'
            'AACTCCATGTGTAGCGGTGGAATGCGTAGATATATGGAAGAACACCAGTGGCGAAGGCGGCTGTC'
            'TAGTCTGTAACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGATACCCTGGTAGT'
            'CCATGCCGTAAACGATGAGTGCTAGGTGTTGGAGGGTTTCCGCCCTTCAGTGCCGCAGCTAACGC'
            'ATTAAGCACTCCGCCTGGGGAGTACGACCGCAAGGTTGAAACTCAAAGGAATTGACGGGGGCCCG'
            'CACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCTACGCGAAGAACCTTACCAGGTCTTGACATA'
            'CTATGCAAACCTAAGAGATTAGGCGTTCCCTTCGGGGACATGGATACAGGTGGTGCATGGTTGTC'
            'GTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATTATCAGTTGC'
            'CAGCATTCAGTTGGGCACTCTGGTGAGACTGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGT'
            'CAAATCATCATGCCCCTTATGACCTGGGCTACACACGTGCTACAATGGACGGTACAACGAGTTGC'
            'GAAGTCGTGAGGCTAAGCTAATCTCTTAAAGCCGTTCTCAGTTCGGATTGTAGGCTGCAACTCGC'
            'CTACATGAAGTTGGAATCGCTAGTAATCGCGGATCAGCATGCCGCGGTGAATACGTTCCCGGGCC'
            'TTGTACACACCGCCCGTCACACCATGAGAGTTTGTAACACCCAAAGCCGGTGAGATAACCTTCGG'
            'GAGTCAGCCGTCTATAGTG',
            '>192680',
            'AGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAACGAAGC'
            'CTTCTTTCACCGAATGTTTGCATTCACCGAAAGAAGCTTAGTGGCGAACGGGTGAGTAACACGTA'
            'GGCAACCTGCCCAAAAGAGGGGGATAACACTTGGAAACAGGTGCTAATACCGCATAACCATGAAC'
            'ACCGCATGATGTTCATGTAAAAGGCGGCTTTTGCTGTCACTTTTGGATGGGCCTGCGGCGTATTA'
            'ACTTGTTGGTGGGGTAACGGCCTACCAAGGTGATGATACGTAGCCGAACTGAGAGGTTGATCGGC'
            'CACATTGGGACTGAGACACGGCCCAAACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATG'
            'GACGAAAGTCTGATGGAGCAACGCCGCGTGAATGAAGAAGGCCTTCGGGTCGTAAAATTCTGTTG'
            'TCAGAGAAGAACGTGCGTGAGAGTAACTGTTCACGTATTGACGGTATCTGATCAGAAAGCCACGG'
            'CTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCGAGCGTTGTCCGGATTTATTGGGCGT'
            'AAAGGGAACGCAGGCGGTCTTTTAAGTCTGATGTGAAAGCCTTCGGCTTAACCGAAGTAGTGCAT'
            'TGGAAACTGGAAGACTTGAGTGCAGAAGAGGAGAGTGGAACTCCATGTGTAGCGGTGAAATGCGT'
            'AGATATATGGAAGAACACCAGTGGCGAAAGCGGCTCTCTGGTCTGTAACTGACGCTGAGGTTCGA'
            'AAGCGTGGGTAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTGCTAAGT'
            'GTTGGAGGGTTTCCGCCCTTCAGTGCTGCAGCTAACGCATTAAGCACTCCGCCTGGGGAGTACGG'
            'TCGCAAGACTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAAT'
            'TCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCTTCTGACAATTCTAGAGATGGAACGTT'
            'CCCTTCGGGGACAGAATGACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTCGTGAGATGTTGGGT'
            'TAAGTCCCGCAACGAGCGCAACCCTTATTGTCAGTTGCCATCATTAAGTTGGGCACTCTGGCGAG'
            'ACTGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAATCATCATGCCCCTTATGACCTGG'
            'GCTACACACGTGCTACAATGGACGGTACAACGAGTCGCTAACTCGCGAGGGCAAGCTAATCTCTT'
            'AAAGCCGTTCTCAGTTCGGACTGCAGGCTGCAACCCGCCTGCACGAAGTTGGAATTGCTAGTAAT'
            'CGCGGATCAGCATGCCGCGGTGAATACGTTCCCGGGTCTTGCACTCACCGCCCGTCA',
            '>2562098',
            'AGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGGATC'
            'ATCGGGAGCTTGCTCCCGATGATCAGCGGCGGACGGGTGAGTAACACGTGGGCAACCTGCCTGTA'
            'AGACTGGGATAACTCCGGGAAACCGGGGCTAATACCGGATAATTCATCTCCTCTCATGAGGGGAT'
            'GCTGAAAGACGGTTTCGGCTGTCACTTACAGATGGGCCCGCGGCGCATTAGCTAGTTGGTGAGGT'
            'AACGGCTCACCAAGGCAACGATGCGTAGCCGACCTGAGAGGGTGATCGGCCACACTGGGACTGAG'
            'ACACGGCCCAGACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCGCAATGGACGAAAGTCTGACG'
            'GAGCAACGCCGCGTGAGCGAAGAAGGCCTTCGGGTCGTAAAGCTCTGTTGTCAGGGAAGAACAAG'
            'TACCGGAGTAACTGCCGGTACCTTGACGGTACCTGACCAGAAAGCCACGGCTAACTACGTGCCAG'
            'CAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGC'
            'GGTCCTTTAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGTCATTGGAAACTGGGGGAC'
            'TTGAGTGCAGAAGAGGAGAGCGGAATTCCACGTGTAGCGGTGAAATGCGTAGAGATGTGGGGGAA'
            'CACCAGTGGCGAAGGCGGCTCTCTGGTCTGTAACTGACGCTGAGGCGCGAAAGCGTGGGGAGCGA'
            'ACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTGCTAAGTGTTAGAGGGTTTCCG'
            'CCCTTTAGTGCTGCAGCAAACGCATTAAGCACTCCGCCTGGGGAGTACGGCCGCAAGGCTGAAAC'
            'TCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAA'
            'GAACCTTACCAGGTCTTGACATCCTCTGCCACTCCTGGAGACAGGACGTTCCCCTTCGGGGGACA'
            'GAGTGACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAAC'
            'GAGCGCAACCCTTGTTCTTAGTTGCCAGCATTCAGTTGGGCGCTCTAAGGAGACTGCCGGTGACA'
            'AACCGGAGGAAGGTGGGGATGACGTCAAATCATCATGCCCCTTATGACCTGGGCTACACACGTGC'
            'TGCAATGGATAGAACAAAGGGCAGCGAAGCCGCGAGGTGAAGCCAATCCCATAAATCTATTCTCA'
            'GTTCGGATTGCAGGCTGCAACTCGCCTGCATGAAGCCGGAATCGCTAGTAATCGCGGATCAGCAT'
            'GCCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCACGAGAGTTTGTAACAC'
            'CCGAAGTCGGTGGGGTAACCTTTTGGAGCCAGCCGCCTAAGGTGGGACAGATGATTGGGGTGAAG'
            'TCGTAACAAGGTA',
            '>4308624',
            'CTTTATTGGAGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAATACATGCAAGTC'
            'GAGCGAATGGATTGAGAGCTTGCTCTCAAGAAGTTAGCGGCGGACGGGTGAGTAACACGTGGGTA'
            'ACCTGCCCATAAGACTGGGATAACTCCGGGAAACCGGGGCTAATACCGGATAACATTTTGAACCG'
            'CATGGTTCGAAATTGAAAGGCGGCTTTGGCTGTCACTTATGGATGGACCCGCGTCGCATTAGCTA'
            'GTTGGTGAGGTAACGGCTCACCAAGGCAACGATGCGTAGCCGACCTGAGAGGGTGATCGGCCACA'
            'CTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCGCAATGGACG'
            'AAAGTCTGACGGAGCAACGCCGCGTGAGTGATGAAGGCTTTCGGGTCGTAAAACTCTGTTGTTAG'
            'GGAAGAACAAGTGCTAGTTGAATAAGCTGGCACCTTGACGGTACCTAACCAGAAAGCCACGGCTA'
            'ACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTATCCGGAATTATTGGGCGTAAA'
            'GCGCGCGCAGGTGGTTTCTTAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGTCATTGG'
            'AAACTGGGAGACTTGAGTGCAGAAGAGGAAAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGA'
            'GATATGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGTAACTGACACTGAGGCGCGAAAG'
            'CGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTGCTAAGTGTT'
            'AGAGGGTTTCCGCCCTTTAGTGCTGAAGTTAACGCATTAAGCACTCCGCCTGGGGAGTACGGCCG'
            'CAAGGCTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCG'
            'AAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCTCTGAAAACCCTAGAGATAGGGCTTCTCC'
            'TTCGGGAGCAGAGTGACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAA'
            'GTCCCGCAACGAGCGCAACCCTTGATCTTAGTTGCCATCATTAAGTTGGGCACTCTAAGGTGACT'
            'GCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAATCATCATGCCCCTTATGACCTGGGCT'
            'ACACACGTGCTACAATGGACGGTACAAAGAGCTGCAAGACCGCGAGGTGGAGCTAATCTCATAAA'
            'ACCGTTCTCAGTTCGGATTGTAGGCTGCAACTCGCCTACATGAAGCTGGAATCGCTAGTAATCGC'
            'GGATCAGCATGCCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCACGAGAG'
            'TTTGTAACACCCGAAGTCGGTGGGGTAACCTTTATGGAGCCAGCCGCCTAAGGTGGGACAGATGA'
            'TTGGGGTGAAGTCGTAACAAGGTAGCCGTATCGGAAGGTGCGGCTGGATCACCTCCTTTCT',
            '>102222',
            'GCGAACGGGTGAGTAACAGGTGGGTACCTGCCCAGAAGCAGGGGATAACACTTGGAAACAGATGT'
            'TAATCCCGTATAACAAAAGAAACCCGCTTGTTTTTCTTTAAAAAGATGGTTGTGCTTATCACTTT'
            'TGATGGACCCGGGGCGCATTAGCTAGTTGGTGAGGTAACGGCTCACCAAGGCAATGATGCGTAGC'
            'CGACTTGAGAGGGTAATCGCCACATTGGGATTGAGACACGGCCCAGACTCTACGGGAGCAGCAGT'
            'AGGGAATCTTCCACCAATGGACGCAAGTCTGATGGAGCAACGCCGCGTGAGTGAAGAAGGGTTTC'
            'GGCTCGTAAAGCTCTGTTGTTAAAGAAGAACGTGGGTGAGAGTAACTGTTCACCCAGTGACGGTA'
            'TTTAACCAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTT'
            'ATCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTCTTTTAAGTCTAATGTGAAAGCCTTCGG'
            'CTCAACCGAAGAAGTGCATTGGAAACTGGGAGACTTGAGTGCAGAAGAGGACAGTGGAACTCCAT'
            'GTGTAGCGGTGAAATGCGTAGATATATGGAAGAACACCAGTGGCGAAGGCGGCTGTCTGGTCTGT'
            'AACTGACGCTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGATACCCTGGTAGTCCATGCCG'
            'TAAACGATGATTACTAAGTGTTGGAGGGTTTCCGCCCTTCAGTGCTGCAGCTAACGCATTAAGTA'
            'ATCCGCCTGGGGAGTACGACCGCAAGGTTGAAACTCAAAAGAATTGACGGGGGCCCGCACAAGCG'
            'GTGGAGCATGTGGTTTAATTCGAAGCTACGCGAAGAACCTTACCAGGTCTTGACATCTTCTGCCA'
            'ACCTAAGAGATTAGGCGTTCCCTTCGGGGACAGAATGACAGGTGGTGCATGGTTGTCTTCAGCTC'
            'GTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATTACTAGTTGCCAGCATTC'
            'AGTTGGGCACTCTAATGAGACTGCCGGTGACAAACCGGAGGAAGGTGGGGACGACGTCAAATCAT'
            'CATGCCCCTTATGACCTGGGCTACACACGTGCTACAATGGATGGGCAACGAGTCCCGAAACCGCG'
            'AGGTTAACTAATCTCTTAAAACCATTCTCAATTCGGACTGTAGGCTGCAACTCGCCTACACGAAG'
            'TCGGAATCGCTAGTAATCGCGGATCAACATGCCGCGGGGAATACGTTCCCGGGCCTTGTACACAC'
            'CGCCGTCACACCATGAGAATTTGTACACCCAAAGCCGGTGGGGTACCTTTTAGACTACCGGCTAA'
            'AGTGGGGACAAATGATTAAGGTGAAGTCGTA',
            '>27815',
            'TAACCCGTGGGCACCTTGCCAAGAAGCAGGGGATAACCCTTGGAAAAGGTTGCTAATTCCGTATA'
            'ACAGAGAAAACCGCCTTGGTTTTCTTTTAAAAGGATGGTTCTGCTATCACTTCTGGATGGACCCG'
            'CGGCGCATTAGCTAGTTGGTGAGGTAACGGCTCACCAAGGCGATGATGCGTAGCCGACCTGAGAG'
            'GTAACCGCCACATTGGACTGAGACACGGGCCAAGACCTCTACGGGAGGCAGAAGTAGGGAATCTT'
            'CGACAATGGACGAAAGTCTGATGGAGCAACGCCGCGTGAGTGAAGAAGGGTTTCGGATCGTAAAG'
            'CTCTGTTGTTAAAGAAGAACGTGGGTGAGAGTAACTGTTCACCCATTGACGGTATTTAACCAGAA'
            'AGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAACCGTTGATCCGGGATT'
            'TATTGGGCGTAAAGCGAGCGCAGGCGGTCTTTTAAGTCTAATGTTGAAAACTTCGGCTCAACCGA'
            'AGAAGTGCATTGGAAACTGGGAGACTTGAGTGCAGAAGAGGATAGTGGAACTCCATGTGTAGCGG'
            'TGAAATGCGTAGATATATGGAGGAACACCAGTGGCGAAGGCGGCTGCTCTGGTCTGTAACTGACG'
            'CTGAGGCTCGAAAGCATGGGTAGCGAACAGGATTAGATACCCTGGTAGTCCATGCCGTAAACGAT'
            'GATTACTAAGTGTTGGAGGGTTTCCGCCCTTCAGTGCTGCAGCTAACGCATTAAGTAATCCGCCT'
            'GGGGAGTACGACCGCAAGGTTGAAACTCAAAAGAATTGACGGGGGCCCGCACAAGCGGTGGAGCA'
            'TGTGGTTTAATTCGAAGCTACGCGAAGAACCTTACCAGGTCTTGACATCTTCTGCCAACCTAGAG'
            'ATTAGGCGTTCCCTTCGGGGACAGAATGACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTCGTGA'
            'GATGTTGGGTTAAGTCCCGCAACGGCGCAACCCTTATTACTAGTTGCCAGCATTCAGTTGGGCAC'
            'TCTAGTGAGACTGCCGGTGACAAACCGGAGGAAGGTGGGGACGACGTCAAATCATCATGCCCCTT'
            'ATGACCTGGGCTACACACGTGCTACAATGGGTGGTACAACGAGTTGCGAAACCGCGAGGTTTAAG'
            'CTAATCTCTTAAAACCATTCTCAGTTCGGACTGTAGGCTGCAACTCGCCTACACGAAGTCGGAAT'
            'CGCTAGTAATCGCGGATCAGCATGCCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGGCCGT'
            'CACACCATGAGAGTTTGTAACAACCCAAGGCGGTTGGGTAACCTTTTAGGGGCTAGCCGTTTAAG'
            'GTGGGACAAATTATTAGGGGTGAAGTCGTAACAAGGTAACC',
            '>1136710',
            'GATGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGGACGCACTTTCGTTGATTGAATTAGAGA'
            'TGCTTGCATCGAAGATGATTTCAACTATAAAGTGAGTGGCGAACGGGTGAGTAACACGTGGGTAA'
            'CCTGCCCAGAAGTGGGGGATAACACCTGGAAACAGATGCTAATACCGCATAATAAAATGAACCGC'
            'ATGGTTTATTTTTAAAAGATGGCTTCGGCTATCACTTCTGGATGGACCCGCGGCGTATTAGCTAG'
            'TTGGTGAGATAAAGGCTCACCAAGGCTGTGATACGTAGCCGACCTGAGAGGGTAATCGGCCACAT'
            'TGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGA'
            'AAGTCTGATGGAGCAACGCCGCGTGAGTGATGAAGGCTTTAGGGTCGTAAAACTCTGTTGTTGGA'
            'GAAGAACGTGTGTGAGAGTAACTGCTCATGCAGTGACGGTATCCAACCAGAAAGCCACGGCTAAC'
            'TACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTATCCGGATTTATTGGGCGTAAAGC'
            'GAGCGCAGGCGGTTTTTTAAGTCTAATGTGAAAGCCTTCGGCTTAACCGAAGAAGTGCATTGGAA'
            'ACTGGGAAACTTGAGTGCAGAAGAGGACAGTGGAACTCCATGTGTAGCGGTGAAATGCGTAGATA'
            'TATGGAAGAACACCAGTGGCGAAGGCGGCTGTCTGGTCTGTAACTGACGCTGAGGCTCGAAAGCA'
            'TGGGTAGCGAACAGGATTAGATACCCTGGTAGTCCATGCCGTAAACGATGAATGCTAAGTGTTGG'
            'AGGGTTTCCGCCCTTCAGTGCTGCAGCTAACGCATTAAGCATTCCGCCTGGGGAGTACGACCGCA'
            'AGGTTGAAACTCAAAAGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAA'
            'GCTACGCGAAGAACCTTACCAGGTCTTGACATCTTCTGCTAACCTAAGAGATTAGGCGTTCCCTT'
            'CGGGGACAGAATGACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGT'
            'CCCGCAACGAGCGCAACCCCTATTATTAGTTGCCAGCATTAAGTTGGGCACTCTAGTGAGACTGC'
            'CGGTGACAAACCGGAGGAAGGTGGGGACGACGTCAAATCATCATGCCCCTTATGACCTGGGCTAC'
            'ACACGTGCTACAATGGACGGTACAACGAGTTGCGAGACCGCGAGGTTAAGCTAATCTCTTAAAAC'
            'CGTTCTCAGTTCGGACTGCAGGCTGCAACTCGCCTGCACGAAGTTGGAATCGCTAGTAATCGCGG'
            'ATCAGCATGCCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGAGAGTT'
            'TGTAACACCCAAAGCCGGTGGAGTAACCTTCGGGAGCTAGCCGTCTAAGGTGGGACAGATGATTG'
            'GGGTGAAGTCGTAACAAGGTAGCCGTAGGAGAACCTGCGGCTGGATCACCTCCTTTCT'])

        _taxa1 = '\n'.join([
            '179419	k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales'
            '; f__Lactobacillaceae; g__Lactobacillus; s__brevis',
            '1117026	k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacil'
            'lales; f__Lactobacillaceae; g__Lactobacillus; s__brevis',
            '192680	k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales'
            '; f__Lactobacillaceae; g__Lactobacillus; s__ruminis',
            '2562098	k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales;'
            ' f__Bacillaceae; g__Bacillus; s__foraminis',
            '4308624	k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales;'
            ' f__Bacillaceae; g__Bacillus; s__cereus',
            '102222	k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales'
            '; f__Lactobacillaceae; g__Pediococcus; s__acidilactici',
            '27815	k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales'
            '; f__Lactobacillaceae; g__Pediococcus; s__acidilactici',
            '1136710	k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacil'
            'lales; f__Lactobacillaceae; g__Pediococcus; s__damnosus'])

        _taxa2 = '\n'.join([
            '179419	k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales'
            '; f__Lactobacillaceae; g__Lactobacillus; s__brevis',
            '1117026	k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacil'
            'lales; f__Lactobacillaceae; g__Lactobacillus',
            '192680	k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales'
            '; f__Lactobacillaceae',
            '2562098	k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales',
            '4308624	k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales;'
            ' f__Bacillaceae; g__Bacillus; s__cereus',
            '102222	k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales;'
            ' f__Bacillaceae; g__Bacillus; s__cereus',
            '27815	k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales'
            '; f__Lactobacillaceae; g__Pediococcus',
            '1136710	k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacil'
            'lales; f__Lactobacillaceae; g__Pediococcus'])

        self.tmpdir = mkdtemp()
        name = 'B1-REF-L6-iter0'
        self.query_fp = join(self.tmpdir, name, 'query_taxa.tsv')
        self.paramdir = join(self.tmpdir, name, name, 'method1', 'param1')
        self.obs_taxa_fp = join(self.paramdir, 'query_tax_assignments.txt')
        refpath = join(self.tmpdir, 'ref1.txt')
        self.cvdir = join(self.tmpdir, 'cross-validated')
        self.ntdir = join(self.tmpdir, 'novel-taxa-simulations')
        if not exists(self.paramdir):
            makedirs(self.paramdir)
        with open(refpath, 'w') as out:
            out.write(_ref1)
        with open(self.query_fp, 'w') as out:
            out.write(_taxa1)
        with open(self.obs_taxa_fp, 'w') as out:
            out.write(_taxa2)

        self.databases = {'B1-REF': [refpath, self.query_fp,
                                     "ref1", "GTGCCAGCMGCCGCGGTAA",
                                     "GGACTACHVGGGTWTCTAAT", "515f", "806r"]}
        self.ref_data = pd.DataFrame.from_dict(self.databases, orient="index")
        self.ref_data.columns = ["Reference file path", "Reference tax path",
                                 "Reference id", "Fwd primer", "Rev primer",
                                 "Fwd primer id", "Rev primer id"]

        self.exp_p = [0, 1.0, 1.0, 0.875, 0.8571428571428571,
                      0.83333333333333337, 0.66666666666666663]
        self.exp_r = [0, 1.0, 1.0, 0.875, 0.75, 0.625, 0.25]
        self.exp_f = [0, 1.0, 1.0, 0.875, 0.79999999999999993,
                      0.7142857142857143, 0.36363636363636365]
        self.exp_m = [0, 0, 0, 1, 1, 1, 3, 2]

    @classmethod
    def tearDownClass(self):
        rmtree(self.tmpdir)


if __name__ == "__main__":
    main()
