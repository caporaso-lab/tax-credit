#!/usr/bin/env python
from __future__ import division
from glob import glob
from os.path import abspath, join, exists, split
from collections import defaultdict

from numpy import asarray
from pylab import scatter, xlabel, ylabel, xlim, ylim

from biom.exception import UnknownIDError

from cogent.maths.stats.test import correlation_test
from cogent.maths.distance_transform import dist_bray_curtis
from cogent.draw.distribution_plots import generate_box_plots
from qiime.transform_coordinate_matrices import procrustes_monte_carlo,\
    get_procrustes_results
from qiime.principal_coordinates import pcoa
from qiime.format import format_distance_matrix

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2013, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

"""Contains code to support the automated evaluation framework.
"""

def find_and_process_result_tables(start_dir,
                                   biom_processor=abspath,
                                   filename_pattern='table*biom'):
    """ given a start_dir, return list of tuples describing the table and containing the processed table

         start_dir: top-level directory to use when starting the walk
         biom_processor: takes a relative path to a biom file and does
          something with it. default is call abspath on it to convert the
          relative path to an absolute path, but could also be
          parse_biom_table, for example. Not sure if we'll want this, but
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
                                     filename_pattern='table*biom'):
    """ given a start_dir, return list of tuples describing the table and containing the processed table

         start_dir: top-level directory to use when starting the walk
         biom_processor: takes a relative path to a biom file and does
          something with it. default is call abspath on it to convert the
          relative path to an absolute path, but could also be
          parse_biom_table, for example. Not sure if we'll want this, but
          it's easy to hook up.
        filename_pattern: pattern to use when matching filenames, can contain
         globbable (i.e., bash-style) wildcards (default: "table*biom")

        results = [(data-set-id, reference-id, biom_processor(table_fp)),
                   ...
                  ]
    """
    table_fps = glob(join(start_dir,'*','*','expected',filename_pattern))
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
                               filename_pattern='table*biom'):
    """ given a start_dir, return list of tuples describing the expected table and containing the processed table

         start_dir: top-level directory to use when starting the walk
         biom_processor: takes a relative path to a biom file and does
          something with it. default is call abspath on it to convert the
          relative path to an absolute path, but could also be
          parse_biom_table, for example. Not sure if we'll want this, but
          it's easy to hook up.
        filename_pattern: pattern to use when matching filenames, can contain
         globbable (i.e., bash-style) wildcards (default: "table*biom")
    """

    results = defaultdict(dict)
    expected_tables = find_and_process_expected_tables(start_dir,biom_processor,filename_pattern)
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
        result = ';'.join(md['taxonomy'][:level+1])
        return result
    return f

def compute_prfs(result_tables,
                 expected_table_lookup,
                 taxonomy_level=6):
    """ Compute precision, recall, and f-measure for result_tables at taxonomy_level

        result_tables: 2d list of tables to be compared to expected tables,
         where the data in the inner list is:
          [dataset_id, reference_database_id, method_id,
           parameter_combination_id, table_fp]
        expected_table_lookup: 2d dict of dataset_id, reference_db_id to BIOM
         table filepath, for the expected result tables
        taxonomy_level: level to compute results

    """
    for dataset_id, reference_id, method_id, params, actual_table_fp in result_tables:
        ## parse the expected table (unless taxonomy_level is specified, this should be
        ## collapsed on level 6 taxonomy)
        try:
            expected_table_fp = expected_table_lookup[dataset_id][reference_id]
        except KeyError:
            raise KeyError, "Can't find expected table for (%s, %s)." % (dataset_id, reference_id)

        try:
            expected_table = parse_biom_table(open(expected_table_fp,'U'))
        except ValueError:
            raise ValueError, "Couldn't parse BIOM table: %s" % expected_table_fp

        ## parse the actual table and collapse it at the specified taxonomic level
        try:
            actual_table = parse_biom_table(open(actual_table_fp,'U'))
        except ValueError:
            raise ValueError, "Couldn't parse BIOM table: %s" % actual_table_fp
        collapse_by_taxonomy = get_taxonomy_collapser(taxonomy_level)
        actual_table = actual_table.collapse(collapse_by_taxonomy, axis='observation', min_group_size=1)

        ## compute precision, recall, and f-measure and yeild those values
        try:
            p,r,f = compute_prf(actual_table,
                                expected_table)
        except ZeroDivisionError:
            p, r, f = -1., -1., -1.
        yield (dataset_id, reference_id, method_id, params, p, r, f)

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

def compute_pearson_spearman(result_tables,
                             expected_table_lookup,
                             taxonomy_level=6):
    """ Compute pearson and spearman correlations and non-parameteric p-values for a set of results

        result_tables: 2d list of tables to be compared to expected tables,
         where the data in the inner list is:
          [dataset_id, reference_database_id, method_id,
           parameter_combination_id, table_fp]
        expected_table_lookup: 2d dict of dataset_id, reference_db_id to BIOM
         table filepath, for the expected result tables
        taxonomy_level: level to compute results
    """
    ### Start code copied directly from compute_prfs - some re-factoring for re-use is
    ### in order here.
    for dataset_id, reference_id, method_id, params, actual_table_fp in result_tables:
        ## parse the expected table (unless taxonomy_level is specified, this should be
        ## collapsed on level 6 taxonomy)
        try:
            expected_table_fp = expected_table_lookup[dataset_id][reference_id]
        except KeyError:
            raise KeyError, "Can't find expected table for (%s, %s)." % (dataset_id, reference_id)

        try:
            expected_table = parse_biom_table(open(expected_table_fp,'U'))
        except ValueError:
            raise ValueError, "Couldn't parse BIOM table: %s" % expected_table_fp

        ## parse the actual table and collapse it at the specified taxonomic level
        try:
            actual_table = parse_biom_table(open(actual_table_fp,'U'))
        except ValueError:
            raise ValueError, "Couldn't parse BIOM table: %s" % actual_table_fp
        collapse_by_taxonomy = get_taxonomy_collapser(taxonomy_level)
        actual_table = actual_table.collapse(collapse_by_taxonomy, axis='observation', min_group_size=1)
        ### End code copied directly from compute_prfs.

        ## compute spearman and pearson correlations
        actual_vector, expected_vector = get_actual_and_expected_vectors(actual_table,
                                                                         expected_table)
        (pearson_corr_coeff,
         pearson_parametric_p_val,
         pearson_permuted_corr_coeffs,
         pearson_nonparametric_p_val,
         pearson_ci) = \
         correlation_test(actual_vector, expected_vector, method='pearson')
        (spearman_corr_coeff,
         spearman_parametric_p_val,
         spearman_permuted_corr_coeffs,
         spearman_nonparametric_p_val,
         spearman_ci) = \
         correlation_test(actual_vector, expected_vector, method='spearman')
        yield (dataset_id,
               reference_id,
               method_id,
               params,
               pearson_corr_coeff,
               pearson_nonparametric_p_val,
               spearman_corr_coeff,
               spearman_nonparametric_p_val)

def compute_procrustes(result_tables,
                       expected_pc_lookup,
                       taxonomy_level=6,
                       num_dimensions=3,
                       random_trials=999):
    """ Compute Procrustes M2 and p-values for a set of results

        result_tables: 2d list of tables to be compared to expected tables,
         where the data in the inner list is:
          [dataset_id, reference_database_id, method_id,
           parameter_combination_id, table_fp]
        expected_pc_lookup: 2d dict of dataset_id, reference_db_id to principal
         coordinate matrices, for the expected result coordinate matrices
        taxonomy_level: level to compute results
    """
    ### Start code copied ALMOST* directly from compute_prfs - some re-factoring for re-use is
    ### in order here. *ALMOST refers to changes to parser and variable names since expected
    ### is a pc matrix here.

    for dataset_id, reference_id, method_id, params, actual_table_fp in result_tables:
        ## parse the expected table (unless taxonomy_level is specified, this should be
        ## collapsed on level 6 taxonomy)
        try:
            expected_pc_fp = expected_pc_lookup[dataset_id][reference_id]
        except KeyError:
            raise KeyError, "Can't find expected table for (%s, %s)." % (dataset_id, reference_id)

        ## parse the actual table and collapse it at the specified taxonomic level
        try:
            actual_table = parse_biom_table(open(actual_table_fp,'U'))
        except ValueError:
            raise ValueError, "Couldn't parse BIOM table: %s" % actual_table_fp
        collapse_by_taxonomy = get_taxonomy_collapser(taxonomy_level)
        actual_table = actual_table.collapse(collapse_by_taxonomy, axis='observation', min_group_size=1)
        ### End code copied directly from compute_prfs.

        # Next block of code, how do I hate thee? Let me count the ways...
        # (1) dist_bray_curtis doesn't take a BIOM Table object
        # (2) pcoa takes a qiime-formatted distance matrix as a list of lines
        # (3) pcoa return a qiime-formatted pc matrix
        # (4) procrustes_monte_carlo needs to pass through the pc "file" multiple
        #     times, so we actually *need* those the pcs that get passed in to be
        #     lists of lines
        dm = dist_bray_curtis(asarray([v for v in actual_table.iterSampleData()]))
        formatted_dm = format_distance_matrix(actual_table.SampleIds,dm)
        actual_pc = pcoa(formatted_dm.split('\n')).split('\n')
        expected_pc = list(open(expected_pc_fp,'U'))

        ## run Procrustes analysis with monte carlo simulation
        actual_m_squared, trial_m_squareds, count_better, mc_p_value =\
         procrustes_monte_carlo(expected_pc,
                                actual_pc,
                                trials=random_trials,
                                max_dimensions=num_dimensions,
                                sample_id_map=None,
                                trial_output_dir=None)

        yield (dataset_id,
               reference_id,
               method_id,
               params,
               actual_m_squared,
               mc_p_value)


def generate_pr_scatter_plots(query_prf,
                              subject_prf,
                              query_color="b",
                              subject_color="r",
                              x_label="Precision",
                              y_label="Recall"):
    """ Generate scatter plot of precision versus recall for query and subject results

        query_prf: precision, recall, and f-measure values as returned
         from compute_prfs for query data
        subject_prf: precision, recall, and f-measure values as returned
         from compute_prfs for subject data
        query_color: the color of the query points (defualt: blue)
        subject_color: the color of the subject points (defualt: red)
        x_label: x axis label for the plot (default: "Precision")
        y_label: y axis label for the plot (default: "Recall")

    """

    # Extract the query precisions and recalls and
    # generate a scatter plot
    query_precisions = [e[4] for e in query_prf]
    query_recalls = [e[5] for e in query_prf]
    scatter(query_precisions,
            query_recalls,
            c=query_color)

    # Extract the subject precisions and recalls and
    # generate a scatter plot
    subject_precisions = [e[4] for e in subject_prf]
    subject_recalls = [e[5] for e in subject_prf]
    scatter(subject_precisions,
            subject_recalls,
            c=subject_color)

    xlim(0,1)
    ylim(0,1)
    xlabel(x_label)
    ylabel(y_label)


def generate_prf_box_plots(query_prf,
                           subject_prf,
                           metric,
                           x_label="Method"):
    """ Generate box plots for precision, recall, or f-measure

        query_prf: precision, recall, and f-measure values as returned
         from compute_prfs for query data
        subject_prf: precision, recall, and f-measure values as returned
         from compute_prfs for subject data
        metric: metric to generate plots for (choices are "precision",
         "recall", and "f-measure")
        x_label: x axis label for the plot (default: "Method")
    """
    metric_lookup = {"precision":(4,"Precision"),
                     "p":(4,"Precision"),
                     "recall":(5,"Recall"),
                     "r":(5,"Recall"),
                     "f-measure":(6,"F-measure"),
                     "f":(6,"F-measure")}
    try:
        metric_idx, y_label = metric_lookup[metric.lower()]
    except KeyError:
        available_metric_desc = ", ".join(metric_lookup.keys())
        error_msg = "Unknown metric: %s. Available choices are: %s" % (metric, available_metric_desc)
        raise KeyError, error_msg


    distributions_by_method = defaultdict(list)
    for e in subject_prf:
        distributions_by_method[e[2]].append(e[metric_idx])
    for e in query_prf:
        distributions_by_method[e[2]].append(e[metric_idx])

    x_tick_labels, distributions = zip(*distributions_by_method.items())
    generate_box_plots(distributions,
                       x_tick_labels = x_tick_labels,
                       x_label = x_label,
                       y_label = y_label,
                       y_min = 0.0,
                       y_max = 1.0)

def generate_precision_box_plots(query_prf,
                                 subject_prf,
                                 x_label="Method"):
    """ Generate precision box plots

        query_prf: precision, recall, and f-measure values as returned
         from compute_prfs for query data
        subject_prf: precision, recall, and f-measure values as returned
         from compute_prfs for subject data
        x_label: x axis label for the plot (default: "Method")
    """
    generate_prf_box_plots(query_prf, subject_prf,"Precision",x_label)

def generate_recall_box_plots(query_prf,
                              subject_prf,
                              x_label="Method"):
    """ Generate recall box plots

        query_prf: precision, recall, and f-measure values as returned
         from compute_prfs for query data
        subject_prf: precision, recall, and f-measure values as returned
         from compute_prfs for subject data
        x_label: x axis label for the plot (default: "Method")
    """
    generate_prf_box_plots(query_prf, subject_prf,"Recall",x_label)

def generate_fmeasure_box_plots(query_prf,
                                subject_prf,
                                x_label="Method"):
    """ Generate f-measure box plots

        query_prf: precision, recall, and f-measure values as returned
         from compute_prfs for query data
        subject_prf: precision, recall, and f-measure values as returned
         from compute_prfs for subject data
        x_label: x axis label for the plot (default: "Method")
    """
    generate_prf_box_plots(query_prf, subject_prf,"F-measure",x_label)

def generate_prf_table(query_prf,
                       subject_prf,
                       sort_metric="F-measure",
                       num_rows=50):
    """ Generate table of precision, recall, and f-measure data

        query_prf: precision, recall, and f-measure values as returned
         from compute_prfs for query data
        subject_prf: precision, recall, and f-measure values as returned
         from compute_prfs for subject data
        sort_metric: metric to sort rows on (choices are "precision",
         "recall", and "f-measure")
        num_rows: number of rows to include in table (default: 50)
    """
    metric_lookup = {"precision":(4,"Precision"),
                     "p":(4,"Precision"),
                     "recall":(5,"Recall"),
                     "r":(5,"Recall"),
                     "f-measure":(6,"F-measure"),
                     "f":(6,"F-measure")}
    try:
        sort_metric_idx, _ = metric_lookup[sort_metric.lower()]
    except KeyError:
        available_metric_desc = ", ".join(metric_lookup.keys())
        error_msg = "Unknown metric: %s. Available choices are: %s" % (sort_metric, available_metric_desc)
        raise KeyError, error_msg

    precision_idx, precision_title = metric_lookup['precision']
    recall_idx, recall_title = metric_lookup['recall']
    fmeasure_idx, fmeasure_title = metric_lookup['f-measure']

    all_prf = query_prf + subject_prf
    # sort by the (-1 * sort_metric) to avoid having to
    # reverse the list in a second step
    header_format = "{:^15} |{:^12} |{:^12} |{:^12} |{:^15} |{:^30}"
    print header_format.format("Data set",
                               precision_title,
                               recall_title,
                               fmeasure_title,
                               "Method",
                               "Parameters")
    row_format = "{:<15} |{:>12} |{:>12} |{:>12} |{:<15} |{:<30}"
    all_prf.sort(key=lambda x: -x[sort_metric_idx])
    for e in all_prf[:num_rows]:
        data = [e[0],
                '%1.3f' % e[precision_idx],
                '%1.3f' % e[recall_idx],
                '%1.3f' % e[fmeasure_idx],
                e[2],
                e[3]]
        print row_format.format(*data)

def generate_correlation_box_plots(query_pearson_spearman,
                                   subject_pearson_spearman,
                                   metric,
                                   x_label="Method"):
    """ Generate box plots for correlation coefficient

        query_pearson_spearman: pearson and spearman data as returned
         from compute_pearson_spearman for query data
        subject_pearson_spearman: pearson and spearman data as returned
         from compute_pearson_spearman for subject data
        metric: metric to generate plots for (choices are "pearson",
         "spearman")
        x_label: x axis label for the plot (default: "Method")
    """
    metric_lookup = {'pearson':(4,"r"),
                     'spearman':(6,"rho")}

    try:
        metric_idx, y_label = metric_lookup[metric.lower()]
    except KeyError:
        available_metric_desc = ", ".join(metric_lookup.keys())
        error_msg = "Unknown metric: %s. Available choices are: %s" % (metric, available_metric_desc)
        raise KeyError, error_msg

    distributions_by_method = defaultdict(list)
    for e in query_pearson_spearman:
        distributions_by_method[e[2]].append(e[metric_idx])
    for e in subject_pearson_spearman:
        distributions_by_method[e[2]].append(e[metric_idx])

    x_tick_labels, distributions = zip(*distributions_by_method.items())
    generate_box_plots(distributions,
                       x_tick_labels = x_tick_labels,
                       x_label = x_label,
                       y_label = y_label,
                       y_min = -1.0,
                       y_max = 1.0)

def generate_pearson_box_plots(subject_pearson_spearman,
                               query_pearson_spearman,
                               x_label="Method"):
    """ Generate box plots for pearson correlation coefficient

        query_pearson_spearman: pearson and spearman data as returned
         from compute_pearson_spearman for query data
        subject_pearson_spearman: pearson and spearman data as returned
         from compute_pearson_spearman for subject data
        x_label: x axis label for the plot (default: "Method")
    """
    generate_correlation_box_plots(subject_pearson_spearman,
                                   query_pearson_spearman,
                                   "pearson",
                                   x_label)

def generate_spearman_box_plots(subject_pearson_spearman,
                                query_pearson_spearman,
                                x_label="Method"):
    """ Generate box plots for spearman correlation coefficient

        query_pearson_spearman: pearson and spearman data as returned
         from compute_pearson_spearman for query data
        subject_pearson_spearman: pearson and spearman data as returned
         from compute_pearson_spearman for subject data
        x_label: x axis label for the plot (default: "Method")
    """
    generate_correlation_box_plots(subject_pearson_spearman,
                                   query_pearson_spearman,
                                   "spearman",
                                   x_label)

def generate_pearson_spearman_table(query_pearson_spearman,
                                    subject_pearson_spearman,
                                    sort_metric="Pearson",
                                    num_rows=50):
    """ Generate table of pearson and spearman data

        query_pearson_spearman: pearson and spearman data as returned
         from compute_pearson_spearman for query data
        subject_pearson_spearman: pearson and spearman data as returned
         from compute_pearson_spearman for subject data
        sort_metric: metric to sort rows on (choices are "pearson"
         and "spearman")
        num_rows: number of rows to include in table (default: 50)
    """
    metric_lookup = {'pearson':(4,"r"),
                     'spearman':(6,"rho")}

    try:
        sort_metric_idx, _ = metric_lookup[sort_metric.lower()]
    except KeyError:
        available_metric_desc = ", ".join(metric_lookup.keys())
        error_msg = "Unknown metric: %s. Available choices are: %s" % (sort_metric, available_metric_desc)
        raise KeyError, error_msg

    all_pearson_spearman = query_pearson_spearman + subject_pearson_spearman

    pearson_idx, pearson_title = metric_lookup['pearson']
    spearman_idx, spearman_title = metric_lookup['spearman']

    # sort by the (-1 * sort_metric) to avoid having to
    # reverse the list in a second step
    all_pearson_spearman.sort(key=lambda x: -x[sort_metric_idx])

    header_format = "{:^15} |{:^12} |{:^12} |{:^15} |{:^30}"
    print header_format.format("Data set",
                               pearson_title,
                               spearman_title,
                               "Method",
                               "Parameters")
    row_format = "{:<15} |{:>12} |{:>12} |{:<15} |{:<30}"
    for e in all_pearson_spearman[:num_rows]:
        data = [e[0],
                '%1.3f' % e[pearson_idx],
                '%1.3f' % e[spearman_idx],
                e[2],
                e[3]]
        print row_format.format(*data)

def generate_procrustes_table(query_procrustes,
                              subject_procrustes,
                              sort_metric="Procrustes",
                              num_rows=50):
    """ Generate table of pearson and spearman data

        query_procrustes: pearson and spearman data as returned
         from compute_procrustes for query data
        subject_procrustes: pearson and spearman data as returned
         from compute_procrustes for subject data
        sort_metric: metric to sort rows on (choices are "procrustes")
         there is currently only one choice, but leaving this in place to
         maintain consistent interface with related functions
        num_rows: number of rows to include in table (default: 50)
    """
    metric_lookup = {'procrustes':(4,"M^2")}

    try:
        sort_metric_idx, _ = metric_lookup[sort_metric.lower()]
    except KeyError:
        available_metric_desc = ", ".join(metric_lookup.keys())
        error_msg = "Unknown metric: %s. Available choices are: %s" % (sort_metric, available_metric_desc)
        raise KeyError, error_msg

    all_procrustes = query_procrustes + subject_procrustes

    procrustes_idx, procrustes_title = metric_lookup['procrustes']

    # sort by the (-1 * sort_metric) to avoid having to
    # reverse the list in a second step
    all_procrustes.sort(key=lambda x: x[sort_metric_idx])

    header_format = "{:^15} |{:^12} |{:^15} |{:^30}"
    print header_format.format("Data set",
                               procrustes_title,
                               "Method",
                               "Parameters")
    row_format = "{:<15} |{:>12} |{:<15} |{:<30}"
    for e in all_procrustes[:num_rows]:
        data = [e[0],
                '%1.3f' % e[procrustes_idx],
                e[2],
                e[3]]
        print row_format.format(*data)
