#!/usr/bin/env python
from __future__ import division
from glob import glob
from os.path import abspath, join, exists, split
from collections import defaultdict
from biom.parse import parse_biom_table

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
                                   biom_processor=abspath):
    """ given a start_dir, return list of tuples describing the table and containing the processed table
    
         start_dir: top-level directory to use when starting the walk
         biom_processor: takes a relative path to a biom file and does
          something with it. default is call abspath on it to convert the
          relative path to an absolute path, but could also be 
          parse_biom_table, for example. Not sure if we'll want this, but 
          it's easy to hook up.
        
        results = [(data-set-id, reference-id, method-id, parameters-id, biom_processor(table_fp)),
                   ...
                  ]
    """
    
    table_fps = glob(join(start_dir,'*','*','*','*','table.biom'))
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
                                     biom_processor=abspath):
    """ given a start_dir, return list of tuples describing the table and containing the processed table
    
         start_dir: top-level directory to use when starting the walk
         biom_processor: takes a relative path to a biom file and does
          something with it. default is call abspath on it to convert the
          relative path to an absolute path, but could also be 
          parse_biom_table, for example. Not sure if we'll want this, but 
          it's easy to hook up.
        
        results = [(data-set-id, reference-id, biom_processor(table_fp)),
                   ...
                  ]
    """
    table_fps = glob(join(start_dir,'*','*','expected','table.L6-taxa.biom'))
    results = []
    for table_fp in table_fps:
        expected_dir, _ = split(table_fp)
        reference_dir, _ = split(expected_dir)
        dataset_dir, reference_id = split(reference_dir)
        _, dataset_id = split(dataset_dir)
        results.append((dataset_id, reference_id, biom_processor(table_fp)))
    return results

def get_expected_tables_lookup(start_dir,
                               biom_processor=abspath):
    results = defaultdict(dict)
    expected_tables = find_and_process_expected_tables(start_dir,biom_processor)
    for dataset_id, reference_id, processed_table in expected_tables:
        results[dataset_id][reference_id] = processed_table
    return results

def get_observed_observation_ids(table,sample_id=None):
    """ Return the set of observation ids with count > 0 in sample_id
    """
    if sample_id is None:
        sample_idx = 0
    else:
        sample_idx = table.SampleIds.index(sample_id)
        
    def f(v,id_,md):
        return v[sample_idx] > 0.0
    
    filtered_table = table.filterObservations(f=f)
    
    result = set(filtered_table.ObservationIds)
    return result


def compute_prf(actual_table,
                expected_table,
                actual_sample_id=None,
                expected_sample_id=None):
    """ Compute precision, recall, and f-measure based on presence/absence of observations
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

def collapse_by_L6_taxonomy(md):
    result = ';'.join(md['taxonomy'][:6])
    return result

def get_taxonomy_collapser(level):
    def f(md):
        result = ';'.join(md['taxonomy'][:level])
        return result
    return f

def compute_prfs(result_tables,
                 expected_table_lookup,
                 taxonomy_level=6):
    """ Compute p, r, and f for a set of results
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
        actual_table = actual_table.collapseObservationsByMetadata(collapse_by_taxonomy)
        
        ## compute precision, recall, and f-measure and yeild those values
        try:
            p,r,f = compute_prf(actual_table,
                                expected_table)
        except ZeroDivisionError:
            p, r, f = -1., -1., -1.
        yield (dataset_id, reference_id, method_id, params, p, r, f)
        
