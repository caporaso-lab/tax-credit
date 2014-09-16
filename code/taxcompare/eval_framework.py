#!/usr/bin/env python
from __future__ import division

# ----------------------------------------------------------------------------
# Copyright (c) 2014--, taxcompare development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from glob import glob
from os.path import abspath, join, exists, split
from collections import defaultdict

from biom.exception import UnknownIDError, TableException
from biom import load_table
from numpy import asarray, zeros
from pylab import scatter, xlabel, ylabel, xlim, ylim
from scipy.spatial.distance import pdist
from scipy.stats import pearsonr, spearmanr
from skbio import DistanceMatrix
from skbio.draw import boxplots
from skbio.stats.distance import mantel
from mpl_toolkits.axes_grid1 import ImageGrid
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def find_and_process_result_tables(start_dir,
                                   biom_processor=abspath,
                                   filename_pattern='table*biom'):
    """ given a start_dir, return list of tuples describing the table and containing the processed table

         start_dir: top-level directory to use when starting the walk
         biom_processor: takes a relative path to a biom file and does
          something with it. default is call abspath on it to convert the
          relative path to an absolute path, but could also be
          load_table, for example. Not sure if we'll want this, but
          it's easy to hook up.
        filename_pattern: pattern to use when matching filenames, can contain
         globbable (i.e., bash-style) wildcards (default: "table*biom")

        results = [(data-set-id, reference-id, method-id, parameters-id, biom_processor(table_fp)),
                   ...
                  ]
    """

    table_fps = glob(join(start_dir,'*','*','*','*',filename_pattern))
    results = []
    for table_fp in table_fps:
        param_dir, _ = split(table_fp)
        method_dir, param_id = split(param_dir)
        reference_dir, method_id = split(method_dir)
        dataset_dir, reference_id = split(reference_dir)
        _, dataset_id = split(dataset_dir)
        results.append((dataset_id, reference_id, method_id, param_id, biom_processor(table_fp)))
    return results

def find_and_process_expected_tables(start_dir,
                                     biom_processor=abspath,
                                     filename_pattern='table.L%d-taxa.biom',
                                     level=6):
    """ given a start_dir, return list of tuples describing the table and containing the processed table

         start_dir: top-level directory to use when starting the walk
         biom_processor: takes a relative path to a biom file and does
          something with it. default is call abspath on it to convert the
          relative path to an absolute path, but could also be
          load_table, for example. Not sure if we'll want this, but
          it's easy to hook up.
        filename_pattern: pattern to use when matching filenames, can contain
         globbable (i.e., bash-style) wildcards (default: "table*biom")

        results = [(data-set-id, reference-id, biom_processor(table_fp)),
                   ...
                  ]
    """
    filename = filename_pattern % level
    table_fps = glob(join(start_dir,'*','*','expected', filename))
    results = []
    for table_fp in table_fps:
        expected_dir, _ = split(table_fp)
        reference_dir, _ = split(expected_dir)
        dataset_dir, reference_id = split(reference_dir)
        _, dataset_id = split(dataset_dir)
        results.append((dataset_id, reference_id, biom_processor(table_fp)))
    return results

def get_expected_tables_lookup(start_dir,
                               biom_processor=abspath,
                               filename_pattern='table.L%d-taxa.biom',
                               level=6):
    """ given a start_dir, return list of tuples describing the expected table and containing the processed table

         start_dir: top-level directory to use when starting the walk
         biom_processor: takes a relative path to a biom file and does
          something with it. default is call abspath on it to convert the
          relative path to an absolute path, but could also be
          load_table, for example. Not sure if we'll want this, but
          it's easy to hook up.
        filename_pattern: pattern to use when matching filenames, can contain
         globbable (i.e., bash-style) wildcards (default: "table*biom")
    """

    results = defaultdict(dict)
    expected_tables = find_and_process_expected_tables(
        start_dir, biom_processor, filename_pattern, level)
    for dataset_id, reference_id, processed_table in expected_tables:
        results[dataset_id][reference_id] = processed_table
    return results

def get_observed_observation_ids(table,sample_id=None):
    """ Return the set of observation ids with count > 0 in sample_id

        table: the biom table object to analyze
        sample_id: the sample_id to test (default is first sample id in table.SampleIds)
    """
    if sample_id is None:
        sample_id = table.sample_ids[0]

    result = []
    for observation_id in table.observation_ids:
        if table.get_value_by_ids(observation_id, sample_id) > 0.0:
            result.append(observation_id)

    return set(result)


def compute_prf(actual_table,
                expected_table,
                actual_sample_id=None,
                expected_sample_id=None):
    """ Compute precision, recall, and f-measure based on presence/absence of observations

        actual_table: table containing results achieved for query
        expected_table: table containing expected results
        actual_sample_id: sample_id to test (default is first sample id in
         actual_table.SampleIds)
        expected_sample_id: sample_id to test (default is first sample id in
         expected_table.SampleIds)
    """
    actual_obs_ids = get_observed_observation_ids(actual_table,
                                                  actual_sample_id)
    expected_obs_ids = get_observed_observation_ids(expected_table,
                                                    expected_sample_id)

    tp = len(actual_obs_ids & expected_obs_ids)
    fp = len(actual_obs_ids - expected_obs_ids)
    fn = len(expected_obs_ids - actual_obs_ids)

    p = tp / (tp + fp)
    r = tp / (tp + fn)
    f = (2 * p * r) / (p + r)

    return p, r, f

def get_taxonomy_collapser(level):
    """ Returns fn to pass to table.collapse

        level: the level to collapse on in the "taxonomy" observation
         metdata category

    """
    def f(id_, md):
        levels = [l.strip() for l in md['taxonomy'].split(';')]
        result = ';'.join(levels[:level])
        return result
    return f

def filter_table(table, min_count, taxonomy_level):
    def f(data_vector, id_, metadata):
        enough_levels = len(metadata['taxonomy']) >= taxonomy_level
        sufficient_count = data_vector.sum() >= min_count
        return enough_levels and sufficient_count
    return table.filter(f, axis='observation', inplace=False)

def compute_mock_results(result_tables, expected_table_lookup,
                         taxonomy_level=6, min_count=10):
    """ Compute precision, recall, and f-measure for result_tables at taxonomy_level

        result_tables: 2d list of tables to be compared to expected tables,
         where the data in the inner list is:
          [dataset_id, reference_database_id, method_id,
           parameter_combination_id, table_fp]
        expected_table_lookup: 2d dict of dataset_id, reference_db_id to BIOM
         table filepath, for the expected result tables
        taxonomy_level: level to compute results

    """
    results = []
    for dataset_id, reference_id, method_id, params, actual_table_fp in result_tables:
        ## parse the expected table (unless taxonomy_level is specified, this should be
        ## collapsed on level 6 taxonomy)
        try:
            expected_table_fp = expected_table_lookup[dataset_id][reference_id]
        except KeyError:
            raise KeyError, "Can't find expected table for (%s, %s)." % (dataset_id, reference_id)

        try:
            expected_table = load_table(expected_table_fp)
        except ValueError:
            raise ValueError, "Couldn't parse BIOM table: %s" % expected_table_fp

        ## parse the actual table and collapse it at the specified taxonomic level
        try:
            actual_table = load_table(actual_table_fp)
        except ValueError:
            raise ValueError, "Couldn't parse BIOM table: %s" % actual_table_fp

        try:
            actual_table = filter_table(actual_table, min_count, taxonomy_level)
        except TableException:
            # if all data is filtered out, move on to the next table
            continue
        collapse_by_taxonomy = get_taxonomy_collapser(taxonomy_level)
        actual_table = actual_table.collapse(collapse_by_taxonomy, axis='observation', min_group_size=1)

        for sample_id in actual_table.sample_ids:
            ## compute precision, recall, and f-measure
            try:
                p,r,f = compute_prf(actual_table,
                                    expected_table,
                                    actual_sample_id=sample_id,
                                    expected_sample_id=sample_id)
            except ZeroDivisionError:
                p, r, f = -1., -1., -1.

            # compute pearson and spearman
            actual_vector, expected_vector =\
                get_actual_and_expected_vectors(actual_table, expected_table,
                    actual_sample_id=sample_id, expected_sample_id=sample_id)

            pearson_r, pearson_p = pearsonr(actual_vector, expected_vector)
            spearman_r, spearman_p = spearmanr(actual_vector, expected_vector)

            results.append((dataset_id, sample_id, reference_id, method_id, params, p, r, f,
                            pearson_r, pearson_p, spearman_r, spearman_p))

    return pd.DataFrame(results, columns=["Dataset", "SampleID", "Reference", "Method",
                                           "Parameters", "Precision", "Recall",
                                           "F-measure", "Pearson r", "Pearson p",
                                           "Spearman r", "Spearman p"])

def get_actual_and_expected_vectors(actual_table,
                                    expected_table,
                                    actual_sample_id=None,
                                    expected_sample_id=None):
    """ Return vectors of obs counts for obs ids observed in specified samples

        actual_table: table containing results achieved for query
        expected_table: table containing expected results
        actual_sample_id: sample_id to test (default is first sample id in
         actual_table.SampleIds)
        expected_sample_id: sample_id to test (default is first sample id in
         expected_table.SampleIds)
    """
    actual_obs_ids = get_observed_observation_ids(actual_table,
                                                  actual_sample_id)
    expected_obs_ids = get_observed_observation_ids(expected_table,
                                                    expected_sample_id)
    all_obs_ids = list(actual_obs_ids | expected_obs_ids)

    if actual_sample_id is None:
        actual_sample_idx = 0
    else:
        actual_sample_idx = actual_table.index(actual_sample_id, axis="sample")

    if expected_sample_id is None:
        expected_sample_idx = 0
    else:
        expected_sample_idx = expected_table.index(expected_sample_id, axis="sample")

    actual_vector = []
    expected_vector = []
    for obs_id in all_obs_ids:
        try:
            actual_obs_idx = actual_table.index(obs_id, axis="observation")
        except UnknownIDError:
            actual_value = 0.0
        else:
            actual_value = actual_table[actual_obs_idx, actual_sample_idx]
        actual_vector.append(actual_value)

        try:
            expected_obs_idx = expected_table.index(obs_id, axis="observation")
        except UnknownIDError:
            expected_value = 0.0
        else:
            expected_value = expected_table[expected_obs_idx,expected_sample_idx]
        expected_vector.append(expected_value)

    return actual_vector, expected_vector


def distance_matrix_from_table(table, metric='braycurtis'):
    """Compute distances between all pairs of samples in table

        This function was written by Greg Caporaso for scikit-bio. It is
        temporarily here, but is under the BSD license.

        Parameters
        ----------
        table : biom.table.Table
        metric : str
            The name of the scipy pairwise distance (``pdist``) function
            to use when generating pairwise distances.

        Returns
        -------
        skbio.core.distance.DistanceMatrix

        Examples
        --------
        Create a biom Table object containing 10 OTUs and 4 samples. This code was
        pulled from http://biom-format.org/documentation/table_objects.html

        >>> import numpy as np
        >>> from biom.table import Table
        >>> data = np.arange(40).reshape(10, 4)
        >>> data[2,2] = 0
        >>> sample_ids = ['S%d' % i for i in range(4)]
        >>> observ_ids = ['O%d' % i for i in range(10)]
        >>> sample_metadata = [{'environment': 'A'}, {'environment': 'B'},
        ...                    {'environment': 'A'}, {'environment': 'B'}]
        >>> observ_metadata = [{'taxonomy': ['Bacteria', 'Firmicutes']},
        ...                    {'taxonomy': ['Bacteria', 'Firmicutes']},
        ...                    {'taxonomy': ['Bacteria', 'Proteobacteria']},
        ...                    {'taxonomy': ['Bacteria', 'Proteobacteria']},
        ...                    {'taxonomy': ['Bacteria', 'Proteobacteria']},
        ...                    {'taxonomy': ['Bacteria', 'Bacteroidetes']},
        ...                    {'taxonomy': ['Bacteria', 'Bacteroidetes']},
        ...                    {'taxonomy': ['Bacteria', 'Firmicutes']},
        ...                    {'taxonomy': ['Bacteria', 'Firmicutes']},
        ...                    {'taxonomy': ['Bacteria', 'Firmicutes']}]
        >>> table = Table(data, observ_ids, sample_ids, observ_metadata,
        ...               sample_metadata, table_id='Example Table')

        Compute Bray-Curtis distances between all pairs of samples and return a
        DistanceMatrix object

        >>> bc_dm = distance_matrix_from_table(table)
        >>> print bc_dm
        4x4 distance matrix
        IDs:
        S0, S1, S2, S3
        Data:
        [[ 0.          0.02702703  0.05263158  0.07692308]
         [ 0.02702703  0.          0.02564103  0.05      ]
         [ 0.05263158  0.02564103  0.          0.02439024]
         [ 0.07692308  0.05        0.02439024  0.        ]]

        Compute Jaccard distances between all pairs of samples and return a
        DistanceMatrix object. (Need a better example here.)

        >>> j_dm = distance_matrix_from_table(table, "jaccard")
        >>> print j_dm
        4x4 distance matrix
        IDs:
        S0, S1, S2, S3
        Data:
        [[ 0.  1.  1.  1.]
         [ 1.  0.  1.  1.]
         [ 1.  1.  0.  1.]
         [ 1.  1.  1.  0.]]

        Determine if the resulting distance matrices are significantly correlated
        by computing the Mantel correlation between them. (Including the p-value
        won't work for doc testing as it's Monte Carlo-based, so exact matching
        will fail.)

        >>> from skbio.math.stats.distance import mantel
        >>> print mantel(j_dm, bc_dm)
        (nan, nan)

        Compute PCoA for both distance matrices, and then find the Procrustes
        M-squared value.
        >>> bc_pc = PCoA(bc_dm).scores()
        >>> j_pc = PCoA(j_dm).scores()
        >>> print procrustes(bc_pc.site, j_pc.site)[2]
        0.645043903715

        Would be really cool to embed a 3d matplotlib scatter plot in here for
        one of the PC matrices... That could make a really cool demo for SciPy.
        I'm thinking one of these:
        http://matplotlib.org/examples/mplot3d/scatter3d_demo.html

    """
    sample_ids = table.sample_ids
    num_samples = len(sample_ids)
    dm = zeros((num_samples, num_samples))
    for i, sid1 in enumerate(sample_ids):
        v1 = table.data(sid1)
        for j, sid2 in enumerate(sample_ids[:i]):
            v2 = table.data(sid2)
            dm[i, j] = dm[j, i] = pdist([v1, v2], metric)
    return DistanceMatrix(dm, sample_ids)

def compute_mantel(result_tables,
                   taxonomy_level=6,
                   random_trials=999):
    """ Compute mantel r and p-values for a set of results

        result_tables: 2d list of tables to be compared,
         where the data in the inner list is:
          [dataset_id, reference_database_id, method_id,
           parameter_combination_id, table_fp]
        taxonomy_level: level to compute results
        random_trials : number of Monte Carlo trials to run in Mantel test
    """
    collapse_by_taxonomy = get_taxonomy_collapser(taxonomy_level)
    results = []

    for dataset_id, reference_id, method_id, params, actual_table_fp in result_tables:
        ## load the table and collapse it at the specified taxonomic level
        try:
            full_table = load_table(actual_table_fp)
        except ValueError:
            raise ValueError("Couldn't parse BIOM table: %s" % actual_table_fp)
        collapsed_table = full_table.collapse(collapse_by_taxonomy,
                                              axis='observation',
                                              min_group_size=1)

        ## Compute Bray-Curtis distances between samples in the full table and
        ## in the collapsed table, and compare them with Mantel.
        # This is way too compute-intensive because we're computing the actual
        # dm everytime, which doesn't need to happen.
        collapsed_dm = distance_matrix_from_table(collapsed_table)
        full_dm = distance_matrix_from_table(full_table)
        mantel_r, p = mantel(collapsed_dm, full_dm)

        results.append((dataset_id, reference_id, method_id, params, mantel_r, p))

    return pd.DataFrame(results, columns=["Dataset", "Reference", "Method",
                                           "Parameters", "Mantel r", "Mantel p"])

def generate_pr_scatter_plots(query_prf,
                              subject_prf,
                              query_color="b",
                              subject_color="r",
                              x_label="Precision",
                              y_label="Recall"):
    """ Generate scatter plot of precision versus recall for query and subject results

        query_prf : pandas.DataFrame
         Precision, recall, and f-measure values as returned from compute_prfs
         for query data
        subject_prf : pandas.DataFrame
         Precision, recall, and f-measure values as returned from compute_prfs
         for subject data
        query_color : str, optional
         The color of the query points
        subject_color : str, optional
         The color of the subject points
        x_label : str, optional
         x axis label for the plot
        y_label : str, optional
         y axis label for the plot

    """
    # Extract the query precisions and recalls and
    # generate a scatter plot
    query_precisions = query_prf['Precision']
    query_recalls = query_prf['Recall']
    scatter(query_precisions,
            query_recalls,
            c=query_color)

    # Extract the subject precisions and recalls and
    # generate a scatter plot
    subject_precisions = subject_prf['Precision']
    subject_recalls = subject_prf['Recall']
    scatter(subject_precisions,
            subject_recalls,
            c=subject_color)

    xlim(0,1)
    ylim(0,1)
    xlabel(x_label)
    ylabel(y_label)

def boxplot_from_data_frame(df,
                            group_by,
                            metric,
                            y_min = 0.0,
                            y_max = 1.0):
    """Generate boxplot of metric by group
    """

    distributions = []
    x_tick_labels = df[group_by].unique()
    x_tick_labels.sort()
    for e in x_tick_labels:
        distribution = df.ix[df[group_by] == e, metric]
        distributions.append(distribution)

    boxplots(distributions,
             x_tick_labels = x_tick_labels,
             x_label = group_by,
             y_label = metric,
             y_min = y_min,
             y_max = y_max)


def heatmap_from_data_frame(df, metric, vmin=0, vmax=1, cmap='Reds'):
    """Generate heatmap of specified metric by (method, parameter) x dataset

    df: pandas.DataFrame

    metric: str
        metric to plot in the heatmap

    """
    # there has to be a better way to build these tuples
    tuples = zip(df["Method"], df["Parameters"])
    index = pd.MultiIndex.from_tuples(tuples=tuples,
                                      names=["Method", "Parameters"])
    df = df.pivot_table(index=index, columns="Dataset", values=metric)
    df.sort()

    height = len(df.index) * 0.35
    width = len(df.columns) * 1

    # Based on SO answer: http://stackoverflow.com/a/18238680
    fig = plt.figure(figsize=(width, height))
    grid = ImageGrid(fig, 111, nrows_ncols=(1, 1),
                     direction='row', axes_pad=0.05, add_all=True,
                     label_mode='1', share_all=False,
                     cbar_location='right', cbar_mode='single',
                     cbar_size='5%', cbar_pad=0.05)

    ax = grid[0]
    ax.set_title(metric, fontsize=20)
    ax.tick_params(axis='both', direction='out', labelsize=12)
    im = ax.imshow(df.values, interpolation='nearest', cmap=cmap,
                   vmax=vmax, vmin=vmin)
    ax.cax.colorbar(im)
    ax.cax.tick_params(labelsize=12)
    ax.set_xticks(np.arange(df.shape[1]))
    ax.set_xticklabels(df.columns.tolist(), rotation=90)
    ax.set_yticks(np.arange(df.shape[0]))
    ax.set_yticklabels(df.index)

    fig.show()
