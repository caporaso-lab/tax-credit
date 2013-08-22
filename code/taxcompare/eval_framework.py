#!/usr/bin/env python
from __future__ import division
from glob import glob
from os.path import abspath, join, exists, split

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
    table_fps = glob(join(start_dir,'*','*','expected','table.biom'))
    results = []
    for table_fp in table_fps:
        reference_dir, _ = split(table_fp)
        dataset_dir, reference_id = split(reference_dir)
        _, dataset_id = split(dataset_dir)
        results.append((dataset_id, reference_id, biom_processor(table_fp)))
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
    
    result = set(table.ObservationIds)
    return result


def compute_prf(actual_table,
                expected_table,
                actual_sample_id=None,
                expected_sample_id=None):
    """
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


