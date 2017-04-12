#!/usr/bin/env python


# ----------------------------------------------------------------------------
# Copyright (c) 2014--, tax-credit development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from sys import exit
from glob import glob
from os.path import abspath, join, exists, split, dirname
from collections import defaultdict
from functools import partial
from random import shuffle
from biom.exception import UnknownIDError, TableException
from biom import load_table
from biom.cli.util import write_biom_table
import pandas as pd
from tax_credit import framework_functions, taxa_manipulator


def get_sample_to_top_params(df, metric, sample_col='SampleID',
                             method_col='Method', dataset_col='Dataset',
                             ascending=False):
    """ Identify the top-performing methods for a given metric

    Parameters
    ----------
    df: pd.DataFrame
    metric: Column header defining the metric to compare parameter combinations
     with
    sample_col: Column header defining the SampleID
    method_col: Column header defining the method name
    dataset_col: Column header defining the dataset name

    Returns
    -------
    pd.DataFrame
     Rows: Multi-index of (Dataset, SampleID)
     Cols: methods
     Values: list of Parameters that achieve the highest performance for each
      method
    """
    sorted_df = df.sort_values(by=metric, ascending=False)
    result = {}

    for dataset in sorted_df[dataset_col].unique():
        dataset_df = sorted_df[sorted_df[dataset_col] == dataset]
        for sid in dataset_df[sample_col].unique():
            dataset_sid_res = dataset_df[dataset_df[sample_col] == sid]
            current_results = {}
            for method in sorted_df.Method.unique():
                m_res = dataset_sid_res[dataset_sid_res.Method == method]
                mad_metric_value = m_res[metric].mad()
                # if higher values are better, find params within MAD of max
                if ascending is False:
                    max_val = m_res[metric].max()
                    tp = m_res[m_res[metric] >= (max_val - mad_metric_value)]
                # if lower values are better, find params within MAD of min
                else:
                    min_val = m_res[metric].min()
                    tp = m_res[m_res[metric] <= (min_val + mad_metric_value)]
                current_results[method] = list(tp.Parameters)
            result[(dataset, sid)] = current_results
    result = pd.DataFrame(result).T
    return result


def parameter_comparisons(
        df, method, metrics=['Precision', 'Recall', 'F-measure',
                             'Taxon Accuracy Rate', 'Taxon Detection Rate'],
        sample_col='SampleID', method_col='Method',
        dataset_col='Dataset', ascending=None):
    """ Count the number of times each parameter combination achieves the top
    score

    Parameters
    ----------
    df: pd.DataFrame
    method: method of interest
    metrics: metrics to include as headers in the resulting DataFrame

    Returns
    -------
    pd.DataFrame
     Rows: Parameter combination
     Cols: metrics, Mean
     Values: Mean: average value of all other columns in row; metrics: count of
      times a parameter combination achieved the best score for the given
      metric
    """
    result = {}

    # Default to descending sort for all metrics
    if ascending is None:
        ascending = {m: False for m in metrics}

    for metric in metrics:
        df2 = get_sample_to_top_params(df, metric, sample_col=sample_col,
                                       method_col=method_col,
                                       dataset_col=dataset_col,
                                       ascending=ascending[metric])
        current_result = defaultdict(int)
        for optimal_parameters in df2[method]:
            for optimal_parameter in optimal_parameters:
                current_result[optimal_parameter] += 1
        result[metric] = current_result
    result = pd.DataFrame.from_dict(result)
    result.fillna(0, inplace=True)
    result = result.sort_values(by=metrics[-1], ascending=False)
    return result


def find_and_process_result_tables(start_dir,
                                   biom_processor=abspath,
                                   filename_pattern='table*biom'):
    """ given a start_dir, return list of tuples describing the table and
    containing the processed table

         start_dir: top-level directory to use when starting the walk
         biom_processor: takes a relative path to a biom file and does
          something with it. default is call abspath on it to convert the
          relative path to an absolute path, but could also be
          load_table, for example. Not sure if we'll want this, but
          it's easy to hook up.
        filename_pattern: pattern to use when matching filenames, can contain
         globbable (i.e., bash-style) wildcards (default: "table*biom")

        results = [(data-set-id, reference-id, method-id, parameters-id,
        biom_processor(table_fp)),
                   ...
                  ]
    """

    table_fps = glob(join(start_dir, '*', '*', '*', '*', filename_pattern))
    results = []
    for table_fp in table_fps:
        param_dir, _ = split(table_fp)
        method_dir, param_id = split(param_dir)
        reference_dir, method_id = split(method_dir)
        dataset_dir, reference_id = split(reference_dir)
        _, dataset_id = split(dataset_dir)
        results.append((dataset_id, reference_id, method_id, param_id,
                        biom_processor(table_fp)))
    return results


def find_and_process_expected_tables(start_dir,
                                     biom_processor=abspath,
                                     filename_pattern='table.L{0}-taxa.biom',
                                     level=6):
    """ given a start_dir, return list of tuples describing the table and
    containing the processed table

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
    filename = filename_pattern.format(level)
    table_fps = glob(join(start_dir, '*', '*', 'expected', filename))
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
                               filename_pattern='table.L{0}-taxa.biom',
                               level=6):
    """ given a start_dir, return list of tuples describing the expected table
    and containing the processed table

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


def get_observed_observation_ids(table, sample_id=None, ws_strip=False):
    """ Return the set of observation ids with count > 0 in sample_id

        table: the biom table object to analyze
        sample_id: the sample_id to test (default is first sample id in
        table.SampleIds)
    """
    if sample_id is None:
        sample_id = table.ids(axis="sample")[0]

    result = []
    for observation_id in table.ids(axis="observation"):
        if table.get_value_by_ids(observation_id, sample_id) > 0.0:
            # remove all whitespace from observation_id
            if ws_strip is True:
                observation_id = "".join(observation_id.split())
            result.append(observation_id)

    return set(result)


def compute_taxon_accuracy(actual_table, expected_table, actual_sample_id=None,
                           expected_sample_id=None):
    """ Compute taxon accuracy and detection rates based on presence/absence of
    observations

        actual_table: table containing results achieved for query
        expected_table: table containing expected results
        actual_sample_id: sample_id to test (default is first sample id in
         actual_table.SampleIds)
        expected_sample_id: sample_id to test (default is first sample id in
         expected_table.SampleIds)
    """
    actual_obs_ids = get_observed_observation_ids(actual_table,
                                                  actual_sample_id,
                                                  ws_strip=True)
    expected_obs_ids = get_observed_observation_ids(expected_table,
                                                    expected_sample_id,
                                                    ws_strip=True)

    tp = len(actual_obs_ids & expected_obs_ids)
    fp = len(actual_obs_ids - expected_obs_ids)
    fn = len(expected_obs_ids - actual_obs_ids)

    if tp > 0:
        p = tp / (tp + fp)
        r = tp / (tp + fn)
    else:
        p, r = 0, 0

    return p, r


def get_taxonomy_collapser(level, md_key='taxonomy',
                           unassignable_label='unassigned'):
    """ Returns fn to pass to table.collapse

        level: the level to collapse on in the "taxonomy" observation
         metdata category

    """
    def f(id_, md):
        try:
            levels = [l.strip() for l in md[md_key].split(';')]
        except AttributeError:
            try:
                levels = [l.strip() for l in md[md_key]]
            # if no metadata is listed for observation, group as Unassigned
            except TypeError:
                levels = [unassignable_label]
        result = ';'.join(levels[:level+1])
        return result
    return f


def filter_table(table, min_count=0, taxonomy_level=None,
                 taxa_to_keep=None, md_key='taxonomy'):
    try:
        _taxa_to_keep = ';'.join(taxa_to_keep)
    except TypeError:
        _taxa_to_keep = None

    def f(data_vector, id_, metadata):
        # if filtering based on number of taxonomy levels, and this
        # observation has taxonomic information, and
        # there are a sufficient number of taxonomic levels
        # Table filtering here removed taxa that have insufficient levels
        enough_levels = taxonomy_level is None or \
                        (metadata[md_key] is not None and
                         len(metadata[md_key]) >= taxonomy_level+1)
        # if filtering to specific taxa, this OTU is assigned to that taxonomy
        allowed_taxa = _taxa_to_keep is None or \
            id_.startswith(_taxa_to_keep) or \
            (metadata is not None and md_key in metadata and
             ';'.join(metadata[md_key]).startswith(_taxa_to_keep))
        # the count of this observation is at least min_count
        sufficient_count = data_vector.sum() >= min_count
        return sufficient_count and allowed_taxa and enough_levels

    return table.filter(f, axis='observation', inplace=False)


def seek_results(results_dirs):
    '''Iterate over a list of directories to find results files and pass these
    to find_and_process_result_tables
    '''
    # Confirm that mock results exist and process tables of observations
    results = []
    for results_dir in results_dirs:
        assert exists(results_dir), '''Mock community result directory
            does not exist: {0}'''.format(results_dir)
        results += find_and_process_result_tables(results_dir)
    return results


def evaluate_results(results_dirs, expected_results_dir, results_fp, mock_dir,
                     taxonomy_level_range=range(2, 7), min_count=0,
                     taxa_to_keep=None, md_key='taxonomy', new_param_ids=None,
                     subsample=False, filename_pattern='table.L{0}-taxa.biom',
                     size=10, per_seq_precision=False, exclude=['other'],
                     force=False):
    '''Load observed and expected observations from tax-credit, compute
        precision, recall, F-measure, and correlations, and return results
        as dataframe.

        results_dirs: list of directories containing precomputed taxonomy
            assignment results to evaluate. Must be in format:
                results_dirs/<dataset name>/
                    <reference name>/<method>/<parameters>/
        expected_results_dir: directory containing expected composition data in
            the structure:
            expected_results_dir/<dataset name>/<reference name>/expected/
        results_fp: path to output file containing evaluation results summary.
        mock_dir: path
            Directory of mock community directiories containing mock feature
            tables without taxonomy.
        taxonomy_level_range: RANGE of taxonomic levels to evaluate.
        min_count: int
            Minimum abundance threshold for filtering taxonomic features.
        taxa_to_keep: list of taxonomies to retain, others are removed before
            evaluation.
        md_key: metadata key containing taxonomy metadata in observed taxonomy
            biom tables.
        new_param_ids: dictionary of lists of parameters for taxonomy
            classifiers to test. In format:
                {classifier_name: [params1, params2]}
        subsample: bool
            Randomly subsample results for test runs.
        size: int
            Size of subsample to take.
        exclude: list
            taxonomies to explicitly exclude from precision scoring.
        force: bool
            Overwrite pre-existing results_fp?
    '''

    # Define the subdirectories where the query mock community data should be
    results = seek_results(results_dirs)

    if subsample is True:
        shuffle(results)
        results = results[:size]

    # Process tables of expected taxonomy composition
    expected_tables = get_expected_tables_lookup(
        expected_results_dir, filename_pattern=filename_pattern)

    # Compute accuracy results OR read in pre-existing mock_results_fp
    if not exists(results_fp) or force is True:
        mock_results = compute_mock_results(results,
                                            expected_tables,
                                            results_fp,
                                            mock_dir,
                                            taxonomy_level_range,
                                            min_count=min_count,
                                            taxa_to_keep=taxa_to_keep,
                                            md_key=md_key,
                                            new_param_ids=new_param_ids,
                                            per_seq_precision=per_seq_precision,  # noqa
                                            exclude=exclude)
    else:
        print("{0} already exists.".format(results_fp))
        print("Reading in pre-computed evaluation results.")
        print("To overwrite, set force=True")
        mock_results = pd.DataFrame.from_csv(results_fp, sep='\t')

    return mock_results


def mount_observations(table_fp, min_count=0, taxonomy_level=6,
                       taxa_to_keep=None, md_key='taxonomy', normalize=True,
                       clean_obs_ids=True, filter_obs=True):
    '''load biom table, filter by abundance, collapse taxonomy, return biom.

    table_fp: path
        Input biom table.
    min_count: int
        Minimum abundance threshold; features detected at lower abundance are
        removed from table.
    taxonomy_level: int
        Taxonomic level at which to collapse table.
    taxa_to_keep: list of taxonomies to retain, others are removed before
        evaluation.
    md_key: str
        biom observation metadata key on which to collapse and filter.
    normalize: bool
        Normalize table to relative abundance across sample rows?
    clean_obs_ids: bool
        Remove '[]()' characters from observation ids? (these are removed from
        the ref db during filtering/cleaning steps, and should be removed from
        expected taxonomy files to avoid mismatches).
    filter_obs: bool
        Filter observations? filter_table will remove observations if taxonomy
        strings are shorter than taxonomic_level, count is less than min_count,
        or observation is not included in taxa_to_keep.
    '''

    try:
        table = load_table(table_fp)
    except ValueError:
        raise ValueError("Couldn't parse BIOM table: {0}".format(table_fp))

    if filter_obs is True and min_count > 0 and taxa_to_keep is not None:
        try:
            table = filter_table(table, min_count, taxonomy_level,
                                 taxa_to_keep, md_key=md_key)
        except TableException:
            # if all data is filtered out, move on to the next table
            pass

        except TypeError:
            print("Missing taxonomic information in table " + table_fp)

        if table.is_empty():
            raise ValueError("Table is empty after filtering at"
                             " {0}".format(table_fp))

    collapse_taxonomy = get_taxonomy_collapser(taxonomy_level, md_key=md_key)

    try:
        table = table.collapse(collapse_taxonomy, axis='observation',
                               min_group_size=1)
    except TableException:
        raise TableException("Failure to collapse taxonomy for table at:"
                             " {0}".format(table_fp))
    except TypeError:
        raise TypeError("Failure to collapse taxonomy: {0}".format(table_fp))

    if normalize is True:
        table.norm(axis='sample')

    return table


def compute_mock_results(result_tables, expected_table_lookup, results_fp,
                         mock_dir, taxonomy_level_range=range(2, 7),
                         min_count=0,
                         taxa_to_keep=None, md_key='taxonomy',
                         new_param_ids=None, per_seq_precision=False,
                         exclude=None):
    """ Compute precision, recall, and f-measure for result_tables at
    taxonomy_level

        result_tables: 2d list of tables to be compared to expected tables,
         where the data in the inner list is:
          [dataset_id, reference_database_id, method_id,
           parameter_combination_id, table_fp]
        expected_table_lookup: 2d dict of dataset_id, reference_db_id to BIOM
         table filepath, for the expected result tables
        taxonomy_level_range: range of levels to compute results
        results_fp: path to output file containing evaluation results summary
        mock_dir: path
            Directory of mock community directories that contain feature tables
            without taxonomy.
        per_seq_precision: bool
            Compute per-sequence precision/recall scores from expected
            taxonomy assignments?
        exclude: list
            taxonomies to explicitly exclude from precision scoring.
    """
    results = []
    new_param_ids = {}
    param_data = {}
    param_ids = {'rdp': ['confidence'], 'blast': ['e-value'],
                 'sortmerna': ['min consensus fraction', 'similarity',
                               'best N alignments', 'coverage', 'e value'],
                 'sortmerna-w16': ['min consensus fraction', 'similarity',
                                   'best N alignments', 'coverage', 'e value'],
                 'uclust': ['min consensus fraction', 'similarity',
                            'max accepts'],
                 'vsearch': ['max accepts', 'similarity',
                             'min consensus fraction'],
                 'blast+': ['e-value', 'max accepts', 'similarity',
                            'min consensus fraction'],
                 'q2-nb': ['confidence']}
    param_ids.update(new_param_ids)
    for dataset_id, ref_id, method, params, actual_table_fp in result_tables:

        # Find expected results
        try:
            expected_table_fp = expected_table_lookup[dataset_id][ref_id]
        except KeyError:
            raise KeyError("Can't find expected table for \
                            ({0}, {1}).".format(dataset_id, ref_id))

        for taxonomy_level in taxonomy_level_range:
            # parse the expected table (unless taxonomy_level is specified,
            # this should be collapsed on level 6 taxonomy)
            expected_table = mount_observations(expected_table_fp,
                                                min_count=0,
                                                taxonomy_level=taxonomy_level,
                                                taxa_to_keep=taxa_to_keep,
                                                filter_obs=False)

            # parse the actual table and collapse it at the specified
            # taxonomic level
            actual_table = mount_observations(actual_table_fp,
                                              min_count=min_count,
                                              taxonomy_level=taxonomy_level,
                                              taxa_to_keep=taxa_to_keep,
                                              md_key=md_key)

            # load the feature table without taxonomy assignment
            # we use this for per-sequence precision
            feature_table_fp = join(mock_dir, dataset_id, 'feature_table.biom')
            try:
                feature_table = load_table(feature_table_fp)
            except ValueError:
                raise ValueError(
                    "Couldn't parse BIOM table: {0}".format(feature_table_fp))

            for sample_id in actual_table.ids(axis="sample"):
                # compute precision, recall, and f-measure
                try:
                    accuracy, detection = compute_taxon_accuracy(
                        actual_table, expected_table,
                        actual_sample_id=sample_id,
                        expected_sample_id=sample_id)
                except ZeroDivisionError:
                    accuracy, detection = -1., -1., -1.

                # compute per-sequence precion / recall
                if per_seq_precision and exists(join(
                        dirname(expected_table_fp), 'trueish-taxonomies.tsv')):
                    p, r, f = per_sequence_precision(
                        expected_table_fp, actual_table_fp, feature_table,
                        sample_id, taxonomy_level, exclude=exclude)
                else:
                    p, r, f = -1., -1., -1.

                # log results
                results.append((dataset_id, taxonomy_level, sample_id,
                                ref_id, method, params, p, r, f, accuracy,
                                detection))

        # record parameter data
        param_data[(method, params)] = {}
        for k, v in zip(param_ids[method], params.split(':')):
            v_ = []
            for e in v:
                try:
                    v_.append(float(e))
                except ValueError:
                    v_.append(e)
            param_data[(method, params)][k] = v_

    param_df = pd.DataFrame(param_data)
    result = pd.DataFrame(results, columns=["Dataset", "Level", "SampleID",
                                            "Reference", "Method",
                                            "Parameters", "Precision",
                                            "Recall", "F-measure",
                                            "Taxon Accuracy Rate",
                                            "Taxon Detection Rate"])
    result = result.merge(param_df.T, left_on=('Method', 'Parameters'),
                          right_index=True)
    result.to_csv(results_fp, sep='\t')
    return result


def per_sequence_precision(expected_table_fp, actual_table_fp, feature_table,
                           sample_id, taxonomy_level, exclude=None):
    '''Precision/recall on individual representative sequences in a mock
    community.
    '''
    # locate expected and observed taxonomies
    exp_fp = join(dirname(expected_table_fp), 'trueish-taxonomies.tsv')
    if exists(exp_fp):
        obs_dir = dirname(actual_table_fp)
        if exists(join(obs_dir, 'rep_seqs_tax_assignments.txt')):
            obs_fp = join(obs_dir, 'rep_seqs_tax_assignments.txt')
        elif exists(join(obs_dir, 'taxonomy.tsv')):
            obs_fp = join(obs_dir, 'taxonomy.tsv')
        else:
            raise RuntimeError('taxonomy assignments do not exist '
                               'for dataset {0}'.format(obs_dir))
        # compile lists of taxa only if observed in current sample
        exp = observations_to_list(exp_fp, feature_table, sample_id)
        obs = observations_to_list(obs_fp, feature_table, sample_id)
        # compile sample weights (observations per sequence in sample)
        weights = [feature_table.get_value_by_ids(
                   line.split('\t')[0], sample_id) for line in exp]
        # load files and run precision/recall
        exp_taxa, obs_taxa = framework_functions.load_prf(
            obs, exp, level=slice(0, taxonomy_level+1))
        ps, rs, fs = framework_functions.compute_prf(
            exp_taxa, obs_taxa, test_type='mock', level=taxonomy_level,
            sample_weight=weights, exclude=exclude)
    else:
        ps, rs, fs = -1., -1., -1.

    return ps, rs, fs


def observations_to_list(obs_fp, actual_table, sample_id):
    '''extract lines from obs_fp to list if they are observed in a given
    sample in biom table actual_table. Returns list of lines from obs_fp,
    which maps biom observation ids (first value in each line) to taxonomy
    labels in tab-delimited file.
    '''
    obs = taxa_manipulator.import_to_list(obs_fp)
    obs = [line for line in obs if
           actual_table.exists(line.split('\t')[0], "observation") and
           actual_table.get_value_by_ids(line.split('\t')[0], sample_id) != 0]
    return obs


def add_sample_metadata_to_table(table_fp, dataset_id, reference_id,
                                 min_count=0, taxonomy_level=6,
                                 taxa_to_keep=None, md_key='taxonomy',
                                 method='expected', params='expected'):
    '''load biom table and populate with sample metadata, then change sample
    names.
    '''

    table = mount_observations(table_fp, min_count=min_count,
                               taxonomy_level=taxonomy_level,
                               taxa_to_keep=taxa_to_keep, md_key=md_key)
    metadata = {s_id: {'sample_id': s_id,
                       'dataset': dataset_id,
                       'reference': reference_id,
                       'method': method,
                       'params': params}
                for s_id in table.ids(axis='sample')}
    table.add_metadata(metadata, 'sample')
    new_ids = {s_id: '_'.join([method, params, s_id])
               for s_id in table.ids(axis='sample')}
    return table.update_ids(new_ids, axis='sample')


def merge_expected_and_observed_tables(expected_results_dir, results_dirs,
                                       md_key='taxonomy', min_count=0,
                                       taxonomy_level=6, taxa_to_keep=None,
                                       biom_fp='merged_table.biom',
                                       filename_pattern='table.L{0}-taxa.biom',
                                       force=False):
    '''For each dataset in expected_results_dir, merge expected and observed
    taxonomy compositions.
    '''
    # Quick and dirty way to keep merge from running automatically in notebooks
    # when users "run all" cells. This is really just a convenience function
    # that is meant to be called from the tax-credit notebooks and causing
    # force=False to kill the function is the best simple control. The
    # alternative is to work out a way to weed out expected_tables that have a
    # merged biom, and just load that biom instead of overwriting if
    # force=False. Then do the same for result_tables. If any new result_tables
    # exist, perform merge if force=True. The only time force=False should
    # result in a new table is when a new mock community/reference dataset
    # combo is added — so just let users set force=True if that's the case.
    if force is False:
        exit('Skipping merge. Set force=True if you intend to generate new '
             'merged tables.')

    # Find expected tables, add sample metadata
    expected_table_lookup = get_expected_tables_lookup(
        expected_results_dir, filename_pattern=filename_pattern)

    expected_tables = {}
    for dataset_id, expected_dict in expected_table_lookup.items():
        expected_tables[dataset_id] = {}
        for reference_id, expected_table_fp in expected_dict.items():
            if not exists(join(expected_results_dir, dataset_id,
                               reference_id, biom_fp)) or force is True:
                expected_tables[dataset_id][reference_id] = \
                    add_sample_metadata_to_table(expected_table_fp,
                                                 dataset_id=dataset_id,
                                                 reference_id=reference_id,
                                                 min_count=min_count,
                                                 taxonomy_level=taxonomy_level,
                                                 taxa_to_keep=taxa_to_keep,
                                                 md_key='taxonomy',
                                                 method='expected',
                                                 params='expected')

    # Find observed results tables, add sample metadata
    result_tables = seek_results(results_dirs)

    for dataset_id, ref_id, method, params, actual_table_fp in result_tables:

        biom_destination = join(expected_results_dir, dataset_id, ref_id,
                                biom_fp)
        if not exists(biom_destination) or force is True:
            try:
                expected_table_fp = \
                    expected_table_lookup[dataset_id][ref_id]
            except KeyError:
                raise KeyError("Can't find expected table for \
                                ({0}, {1}).".format(dataset_id, ref_id))

            # import expected table, amend sample ids
            actual_table = \
                add_sample_metadata_to_table(actual_table_fp,
                                             dataset_id=dataset_id,
                                             reference_id=ref_id,
                                             min_count=min_count,
                                             taxonomy_level=taxonomy_level,
                                             taxa_to_keep=taxa_to_keep,
                                             md_key='taxonomy',
                                             method=method,
                                             params=params)

            # merge expected and resutls tables
            expected_tables[dataset_id][ref_id] = \
                expected_tables[dataset_id][ref_id].merge(actual_table)

            # write biom table to destination
            write_biom_table(expected_tables[dataset_id][ref_id],
                             'hdf5', biom_destination)


def _is_first(df, test_field='Method'):
    """used to filter df to contain only one row per method"""
    observed = set()
    result = []
    for e in df[test_field]:
        result.append(e not in observed)
        observed.add(e)
    return result


def method_by_dataset(df, dataset, sort_field, display_fields,
                      group_by='Dataset', test_field='Method'):
    """ Generate summary of best parameter set for each method for single df
    """
    dataset_df = df.loc[df[group_by] == dataset]
    sorted_dataset_df = dataset_df.sort_values(by=sort_field, ascending=False)
    filtered_dataset_df = sorted_dataset_df[_is_first(sorted_dataset_df,
                                                      test_field)]
    return filtered_dataset_df.ix[:, display_fields]

method_by_dataset_a1 = partial(method_by_dataset,
                               sort_field="F-measure",
                               display_fields=("Method", "Parameters",
                                               "Precision", "Recall",
                                               "F-measure",
                                               "Taxon Accuracy Rate",
                                               "Taxon Detection Rate"))


def method_by_reference_comparison(df, group_by='Reference', dataset='Dataset',
                                   level_range=range(4, 7), lv="Level",
                                   sort_field="F-measure",
                                   display_fields=("Reference", "Level",
                                                   "Method", "Parameters",
                                                   "Precision", "Recall",
                                                   "F-measure",
                                                   "Taxon Accuracy Rate",
                                                   "Taxon Detection Rate")):
    '''Compute mean performance for a given reference/method/parameter
    combination across multiple taxonomic levels.

    df: pandas df
    group_by: str
        Category in df. Means will be averaged across these groups.
    dataset: str
        Category in df. df will be separated by datasets prior to computing
        means.
    level_range: range
        Taxonomy levels to iterate.
    lv: str
        Category in df that contains taxonomic level information.
    sort_field: str
        Category in df. Results within each group/level combination will be
        sorted by this field.
    display_fields: tuple
        Categories in df that should be printed to results table.
    '''

    rank = pd.DataFrame()
    for ds in df[dataset].unique():
        df1 = df[df[dataset] == ds]
        for level in level_range:
            for group in df1[group_by].unique():
                a = method_by_dataset(df1[df1[lv] == level],
                                      group_by=group_by,
                                      dataset=group,
                                      sort_field=sort_field,
                                      display_fields=display_fields)
                rank = pd.concat([rank, a])
    return rank
