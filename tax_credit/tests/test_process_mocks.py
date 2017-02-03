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
from biom import Table
from tax_credit.process_mocks import (amend_biom_taxonomy_ids)


class EvalFrameworkTests(TestCase):

    def setUp(self):
        self.table1 = Table.from_json(json.load(StringIO(table1)))

    def test_amend_biom_taxonomy_ids(self):
        self.assertEqual(set(amend_biom_taxonomy_ids(self.table1,
                             clean_obs_ids=False).ids(axis='observation')),
                         {'k__Archaea;p__;c__;o__;f__;g__;s__',
                         'k__Bacteria;p__;c__;o__;f__;g__;s__',
                         'k__[Fungi];p__;c__;o__;f__;g__;s__'})
        # This also tests clean_taxonomy_ids()
        self.assertEqual(set(amend_biom_taxonomy_ids(self.table1,
                             clean_obs_ids=True).ids(axis='observation')),
                         {'k__Archaea;p__;c__;o__;f__;g__;s__',
                         'k__Bacteria;p__;c__;o__;f__;g__;s__',
                         'k__Fungi;p__;c__;o__;f__;g__;s__'})

table1 = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "greg","date": "2013-08-22T13:10:23.907145","matrix_type": "sparse","matrix_element_type": "float","shape": [3, 4],"data": [[0,0,1.0],[0,1,2.0],[0,2,3.0],[0,3,4.0],[1,0,2.0],[1,1,0.0],[1,2,7.0],[1,3,8.0],[2,0,9.0],[2,1,10.0],[2,2,11.0],[2,3,12.0]],"rows": [{"id": "k__Bacteria", "metadata": {"domain": "Bacteria"}},{"id": "k__Archaea", "metadata": {"domain": "Archaea"}},{"id": "k__[Fungi]", "metadata": {"domain": "[Fungi]"}}],"columns": [{"id": "s1", "metadata": {"country": "Peru", "pH": 4.2}},{"id": "s2", "metadata": {"country": "Peru", "pH": 5.2}},{"id": "s3", "metadata": {"country": "Peru", "pH": 5.0}},{"id": "s4", "metadata": {"country": "Peru", "pH": 4.9}}]}"""

## table 1
# OTU ID	   s1	s2	s3	s4
# k__Archaea    1.0 2.0 3.0 4.0
# k__Bacteria    2.0 0.0 7.0 8.0
# k__[Fungi]    9.0 10.0    11.0    12.0


if __name__ == "__main__":
    main()
