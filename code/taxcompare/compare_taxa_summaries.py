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

from numpy import array, sqrt
from cogent.maths.stats.test import pearson

comparison_modes = ['paired', 'expected']
correlation_types = ['pearson', 'spearman']

def compare_taxa_summaries(taxa_summary1, taxa_summary2, comparison_mode,
                           correlation_type='pearson', normalized_count=1):
    if correlation_type == 'pearson':
        correlation_fn = pearson
    elif correlation_type == 'spearman':
        correlation_fn = spearman_correlation
    else:
        raise ValueError("Invalid correlation type '%s'. Must be one of %r." %
                         (correlation_type, correlation_types))

    taxa_summary1, taxa_summary2 = _sort_and_fill_taxa_summaries(
            [taxa_summary1, taxa_summary2])

    if comparison_mode == 'paired':
        pass
    elif comparison_mode == 'expected':
        correlations = _compare_all_to_expected(taxa_summary1, taxa_summary2,
                                                correlation_fn)
    else:
        raise ValueError("Invalid comparison mode '%s'. Must be one of %r." %
                         (comparison_mode, comparison_modes))

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

def _compare_all_to_expected(observed_taxa_summary, expected_taxa_summary,
                             correlation_fn):
    # Compare each sample in the first taxa summary to the single sample in
    # the second. Make sure that the second taxa summary has only one
    # sample.
    if len(expected_taxa_summary[0]) != 1:
        raise ValueError("The second taxa summary file must contain a single "
                "sample (column) to compare all samples in the first taxa "
                "summary file against when the comparison mode is 'expected'. "
                "You provided %d samples." % len(expected_taxa_summary[0]))

    result = []
    for sample_num in range(len(observed_taxa_summary[0])):
        result.append(correlation_fn(
            observed_taxa_summary[2].T[sample_num],
            expected_taxa_summary[2].T[0]))
#        print "x: %r" % observed_taxa_summary[2].T[sample_num]
#        print "y: %r" % expected_taxa_summary[2].T[0]
    return result

# These next two functions are taken from the qiime.stats.BioEnv class. We
# don't use the class itself because we don't have the required data to
# instantiate one, so the methods are copied here for convenience. TODO: this
# should be moved into pycogent at some point.
def spearman_correlation(vec1, vec2, ranked=False):
    """Calculates the the Spearman distance of two vectors."""
    try:
        temp = len(vec1)
    except ValueError:
        raise ValueError('First input vector is not a list.')
    try:
        temp = len(vec2)
    except ValueError:
        raise ValueError('Second input vector is not a list.')
    if len(vec1) == 0 or len(vec2) == 0:
        raise ValueError('One or both input vectors has/have zero elements')
    if len(vec1) != len(vec2):
        raise ValueError('Vector lengths must be equal')

    if not ranked:
        rank1, ties1 = _get_rank(vec1)
        rank2, ties2 = _get_rank(vec2)
    else:
        rank1, ties1 = vec1
        rank2, ties2 = vec2

    if ties1 == 0 and ties2 == 0:
        n = len(rank1)
        sum_sqr = sum([(x-y)**2 for x,y in zip(rank1,rank2)])
        rho = 1 - (6*sum_sqr/(n*(n**2 - 1)))
        return rho

    avg = lambda x: sum(x)/len(x)

    x_bar = avg(rank1)
    y_bar = avg(rank2)

    numerator = sum([(x-x_bar)*(y-y_bar) for x,y in zip(rank1, rank2)])
    denominator = sqrt(sum([(x-x_bar)**2 for x in rank1])*
                       sum([(y-y_bar)**2 for y in rank2]))
    # Calculate rho.
    return numerator/denominator

def _get_rank(data):
    """Ranks the elements of a list. Used in Spearman correlation."""
    indices = range(len(data))
    ranks = range(1,len(data)+1)
    indices.sort(key=lambda index:data[index])
    ranks.sort(key=lambda index:indices[index-1])
    data_len = len(data)
    i = 0
    ties = 0
    while i < data_len:
        j = i + 1
        val = data[indices[i]]
        try:
            val += 0
        except TypeError:
            raise(TypeError)

        while j < data_len and data[indices[j]] == val:
            j += 1
        dup_ranks = j - i
        val = float(ranks[indices[i]]) + (dup_ranks-1)/2.0
        for k in range(i, i+dup_ranks):
            ranks[indices[k]] = val
        i += dup_ranks
        ties += dup_ranks-1
    return ranks, ties
