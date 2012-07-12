#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout", "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"
__status__ = "Development"

"""Contains functions used in the compare_taxa_summaries.py script."""

from numpy import array
from qiime.parse import parse_taxa_summary_table

def compare_taxa_summaries(taxa_summary_f1, taxa_summary_f2, comparison_mode,
                           correlation_type='pearson', normalized_count=1):
    tax_sum1, tax_sum2 = _sort_and_fill_taxa_summaries([
        parse_taxa_summary_table(taxa_summary_f1),
        parse_taxa_summary_table(taxa_summary_f2)])

    if comparison_mode == 'paired':
        pass
    elif comparison_mode == 'expected':
        pass
    else:
        raise ValueError("Invalid comparison mode '%s'. Must be one of %r." %
                         (comparison_mode, ['paired', 'expected']))

    for i in range(len(taxa_summaries[0][0])):
        for j in range(len(taxa_summaries[1][0])):
            g, p = G_fit((taxa_summary1[i] + pseudo_count) * normalized_count,
                         (taxa_summary2[j] + pseudo_count) * normalized_count)
            yield taxa_summaries[0][0][i], taxa_summaries[1][0][j], g, p

def _sort_and_fill_taxa_summaries(taxa_summaries):
    master_taxa = []
    for ts in taxa_summaries:
        master_taxa += ts[1]
    master_taxa = list(set(master_taxa))
    master_taxa.sort()

    result = []
    for ts in taxa_summaries:
        samples = ts[0]
        orig_taxa = ts[1]
        orig_data = ts[2]
        data = []
        for taxa in master_taxa:
            try:
                taxa_index = orig_taxa.index(taxa)
                data.append(orig_data[taxa_index])
            except ValueError:
                data.append([0.] * len(samples))
        result.append((samples, master_taxa, array(data)))
    return result
