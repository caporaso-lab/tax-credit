#!/usr/bin/env python


# ----------------------------------------------------------------------------
# Copyright (c) 2016--, taxcompare development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os import path, makedirs, remove, rename
from os.path import expandvars, exists, basename, join
from shutil import copyfile
from urllib.request import urlretrieve

from tax_credit.taxa_manipulator import import_taxonomy_to_dict

import qiime
from qiime.plugins import feature_table, demux, dada2, alignment, phylogeny


def extract_mockrobiota_dataset_metadata(mockrobiota_dir, communities):
    '''Extract mock community metadata from mockrobiota dataset-metadata.tsv
    files
    mockrobiota_dir: PATH to mockrobiota directory
    communities: LIST of mock communities to extract
    '''
    dataset_metadata_dict = dict()
    for community in communities:
        dataset_metadata = import_taxonomy_to_dict(join(mockrobiota_dir,
                                                        "data",
                                                        community,
                                                        "dataset-metadata.tsv"
                                                       ))
        dataset_metadata_dict[community] = \
                                (dataset_metadata['raw-data-url-forward-read'],
                                 dataset_metadata['raw-data-url-index-read'],
                                 dataset_metadata['target-gene']
                                )
    return dataset_metadata_dict


def extract_mockrobiota_data(communities, community_metadata, reference_dbs,
                             mockrobiota_dir, mock_data_dir,
                             expected_data_dir):
    '''Extract sample metadata, raw data files, and expected taxonomy

    from mockrobiota, copy to new destination
    communities: LIST of mock communities to extract
    community_metadata: DICT of metadata for mock community.
        see extract_mockrobiota_dataset_metadata()
    reference_dbs = DICT mapping marker_gene to reference set names
    mockrobiota_dir = PATH to mockrobiota repo directory
    mock_data_dir = PATH to destination directory
    expected_data_dir = PATH to destination for expected taxonomy files
    '''
    for community in communities:
        # extract dataset metadata/params
        forward_read_url, index_read_url, marker_gene = \
                                                  community_metadata[community]
        ref_dest_dir, ref_source_dir, ref_version = \
                                                reference_dbs[marker_gene][0:3]

        # mockrobiota source directory
        mockrobiota_community_dir = join(mockrobiota_dir, "data", community)

        # new mock community directory
        community_dir = join(mock_data_dir, community)
        if not exists(community_dir):
            makedirs(community_dir)
        # copy sample-metadata.tsv
        copyfile(join(mockrobiota_community_dir, 'sample-metadata.tsv'),
                 join(community_dir, 'sample-metadata.tsv'))
        # download raw data files
        for file_url in [forward_read_url, index_read_url]:
            destination = join(community_dir, basename(file_url))
            if not exists(destination):
                urlretrieve(file_url, destination)

        # new directory containing expected taxonomy assignments at each level
        expected_taxa_dir = join(expected_data_dir, community,
                                 ref_dest_dir, "expected")
        if not exists(expected_taxa_dir):
            makedirs(expected_taxa_dir)
        # copy expected taxonomy.tsv --- convert to biom???
        copyfile(join(mockrobiota_community_dir, ref_source_dir,
                      ref_version, 'expected-taxonomy.tsv'),
                 join(expected_taxa_dir, 'expected-taxonomy.tsv'))
