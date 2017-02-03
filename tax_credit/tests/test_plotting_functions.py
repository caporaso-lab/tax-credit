#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2014--, tax-credit development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


from unittest import TestCase, main
import json
from io import StringIO
import pandas as pd
from biom.cli.util import write_biom_table
from biom import Table
from os import remove
from tax_credit.plotting_functions import (calculate_linear_regress,
                                           per_level_kruskal_wallis,
                                           make_distance_matrix,
                                           per_method_distance,
                                           within_between_category_distance,
                                           isolate_top_params)


class EvalFrameworkTests(TestCase):

    def setUp(self):
        self.table1 = pd.DataFrame([(n.split('\t'))
                                    for n in table1.split('\n')],
                                   columns=['group', 'dataset','level',
                                            'x', 'y', 'c'],
                                   dtype=float)

    def test_calculate_linear_regress(self):
        result = calculate_linear_regress(self.table1, 'x', 'y', 'group')
        self.assertEqual(result["Slope"][0], 0.5)
        self.assertEqual(result["Slope"][1], 1)
        self.assertEqual(result["Intercept"][0], 0.5)
        self.assertEqual(result["Intercept"][1], 2.0)

    def test_per_level_kruskal_wallis(self):
        result = per_level_kruskal_wallis(
            self.table1, ['c'], 'dataset', 'group', levelrange=range(1, 2))
        self.assertEqual(result[1][0], 0.046301594915110038)

    def test_make_distance_matrix(self):
        self.assertEqual(dm[0][1], 0.16666666666666666)
        self.assertEqual(dm[0][2], 0.27272727272727271)
        self.assertEqual(dm[0][3], 0.33333333333333331)
        self.assertEqual(dm[1][2], 0.27272727272727271)
        self.assertEqual(dm[1][3], 0.33333333333333331)
        self.assertEqual(dm[2][3], 0.066666666666666666)
        self.assertEqual(set(dm.ids), set(s_md.index))

    def test_per_method_distance(self):
        self.assertEqual(dist.iloc[0]['distance'], 0.33333333333333331)
        self.assertEqual(dist.iloc[1]['distance'], 0.33333333333333331)
        self.assertEqual(dist.iloc[2]['distance'], 0.066666666666666666)

    def test_within_between_category_distance(self):
        dist = within_between_category_distance(dm, s_md, 'method',
                                                distance='distance')
        self.assertEqual(dist.iloc[0]['distance'], 0.16666666666666666)
        self.assertEqual(dist.iloc[1]['distance'], 0.27272727272727271)
        self.assertEqual(dist.iloc[2]['distance'], 0.27272727272727271)
        self.assertEqual(dist.iloc[3]['distance'], 0.33333333333333331)
        self.assertEqual(dist.iloc[4]['distance'], 0.33333333333333331)
        self.assertEqual(dist.iloc[5]['distance'], 0.066666666666666666)
        for i in range(3):
            self.assertEqual(dist.iloc[i]['Comparison'], 'within')
            self.assertEqual(dist.iloc[i]['method'], 'A')
        for i in range(3, 6):
            self.assertEqual(dist.iloc[i]['Comparison'], 'between')
            self.assertEqual(dist.iloc[i]['method'], 'B_A')

    def test_isolate_top_params(self):
        best, param_report = isolate_top_params(
            dist, "method", "parameters", "distance", True)
        self.assertEqual(best.iloc[0]['distance'], 0.066666666666666666)
        self.assertEqual(param_report, [('A', 'C')])
        self.assertEqual(len(best), 1)
        best, param_report = isolate_top_params(
            dist, "method", "parameters", "distance", False)
        self.assertEqual(best.iloc[0]['distance'], 0.33333333333333331)

        remove('table2.biom')



table1 = """a\ta\t1\t0.0\t0.5\t0.1
a\ta\t1\t1.0\t1.0\t0.2
a\ta\t1\t2.0\t1.5\t0.2
a\tb\t1\t3.0\t2.0\t8.
a\tb\t1\t4.0\t2.5\t9.
a\tb\t1\t5.0\t3.0\t10.
b\ta\t1\t0.0\t2.0\t0.1
b\ta\t1\t1.0\t3.0\t0.3
b\ta\t1\t2.0\t4.0\t0.1
b\tb\t1\t3.0\t5.0\t9.
b\tb\t1\t4.0\t6.0\t11.
b\tb\t1\t5.0\t7.0\t10."""

table2 = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "greg","date": "2013-08-22T13:10:23.907145","matrix_type": "sparse","matrix_element_type": "float","shape": [3, 4],"data": [[0,0,1.0],[0,1,2.0],[0,2,3.0],[0,3,4.0],[1,0,2.0],[1,1,0.0],[1,2,7.0],[1,3,8.0],[2,0,9.0],[2,1,10.0],[2,2,11.0],[2,3,12.0]],"rows": [{"id": "o1", "metadata": {"domain": "Archaea"}},{"id": "o2", "metadata": {"domain": "Bacteria"}},{"id": "o3", "metadata": {"domain": "Bacteria"}}],"columns": [{"id": "s1", "metadata": {"method": "A", "Sample": "A", "parameters": "A"}},{"id": "s2", "metadata": {"method": "A", "Sample": "A", "parameters": "B"}},{"id": "s3", "metadata": {"method": "A", "Sample": "A", "parameters": "C"}},{"id": "s4", "metadata": {"method": "B", "Sample": "A", "parameters": "D"}}]}"""

## table 2
#OTU ID	s1	s2	s3	s4
# o1    1.0 2.0 3.0 4.0
# o2    2.0 0.0 7.0 8.0
# o3    9.0 10.0    11.0    12.0

table2 = Table.from_json(json.load(StringIO(table2)))
write_biom_table(table2, 'hdf5', 'table2.biom')
dm, s_md = make_distance_matrix('table2.biom', method="braycurtis")
dist = per_method_distance(dm, s_md, group_by='method',
                           standard='B', metric='distance',
                           sample='Sample')


if __name__ == "__main__":
    main()
