#!/usr/bin/env python


# ----------------------------------------------------------------------------
# Copyright (c) 2016--, taxcompare development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os import path, makedirs, remove, rename
from os.path import exists, join
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
        seqs_dir = join(community_dir, 'raw_seqs')
        if not exists(seqs_dir):
            makedirs(seqs_dir)
        # copy sample-metadata.tsv
        copyfile(join(mockrobiota_community_dir, 'sample-metadata.tsv'),
                 join(community_dir, 'sample-metadata.tsv'))
        # download raw data files
        for file_url_dest in [(forward_read_url, 'sequences.fastq.gz'),
                              (index_read_url, 'barcodes.fastq.gz')]:
            destination = join(seqs_dir, file_url_dest[1])
            if not exists(destination):
                urlretrieve(file_url_dest[0], destination)

        # new directory containing expected taxonomy assignments at each level
        expected_taxa_dir = join(expected_data_dir, community,
                                 ref_dest_dir, "expected")
        if not exists(expected_taxa_dir):
            makedirs(expected_taxa_dir)
        # copy expected taxonomy.tsv --- convert to biom???
        copyfile(join(mockrobiota_community_dir, ref_source_dir,
                      ref_version, 'expected-taxonomy.tsv'),
                 join(expected_taxa_dir, 'expected-taxonomy.tsv'))


def batch_demux(communities,
                mock_data_dir,
                demux_params,
                raw_seqs='raw_seqs',
                metadata_fn='sample-metadata.tsv',
                index_col='BarcodeSequence',
                demux_fn='demux-seqs.qza',
                summary_fn='demux_summary.qzv',
                qual_plot_fn='demux_plot_qual.qzv',
                n_qual_plots=1):
    '''raw fastq -> SampleData[SequencesWithQuality]

    Demultiplex raw fastq files, summarize demux, and plot qual scores on
    batch of files

        communities: list
            list of mock communities to extract
        mock_data_dir = filepath
            PATH to data dir containing communities
        demux_params = dict
            DICT of TUPLES listing demux parameters.
            {community : (rev_comp_barcodes, rev_comp_mapping_barcodes)}
        raw_seqs = str
            name of directory containing raw seqs, located in directory
            mock_data_dir/community. Must only contain the following files:
                sequences.fastq.gz
                barcodes.fastq.gz
        metadata_fn = filename
            fn of metadata map, located in mock_data_dir/community
        index_col = str
            column name of metadata column containing index/barcode sequences
        demux_fn = str
            filename to save demux artifact
        summary_fn = str
            filename to save demux summary visualization
        qual_plot_fn = str
            filename to save fastq quality plot visualization
        n_qual_plots = int
            number of fastq quality plots to print
    '''

    for community in communities:
        # extract dataset metadata/params
        community_dir = join(mock_data_dir, community)
        seqs_dir = join(community_dir, raw_seqs)
        sample_md = join(community_dir, metadata_fn)

        # demultiplex
        demux_to_plot_quality(seqs_dir = seqs_dir,
                              sample_md = sample_md,
                              community_dir = community_dir,
                              index_col = index_col,
                              rc_barcodes = demux_params[community][0],
                              rc_map_barcodes = demux_params[community][1],
                              demux_fn = demux_fn,
                              summary_fn = summary_fn,
                              qual_plot_fn = qual_plot_fn,
                              n_qual_plots = n_qual_plots)

        print("{0} complete".format(community))


def demux_to_plot_quality(seqs_dir,
                          sample_md,
                          community_dir,
                          index_col='BarcodeSequence',
                          rc_barcodes=False,
                          rc_map_barcodes=False,
                          demux_fn='demux-seqs.qza',
                          summary_fn='demux_summary.qzv',
                          qual_plot_fn='demux_plot_qual.qzv',
                          n_qual_plots=1):
    '''raw fastq -> SampleData[SequencesWithQuality]

    Demultiplex raw fastq files, summarize demux, and plot qual scores.

        seqs_dir = filepath
            directory containing fastq sequences
        sample_md = filepath
            filepath to sample metadata file
        community_dir: path
            destination directory to print results
        index_col = str
            column name of metadata column containing index/barcode sequences
        rc_barcodes = bool
            reverse complement barcodes prior to demultiplexing?
        rc_map_barcodes = bool
            reverse complement metadata map barcodes prior to demultiplexing?
        demux_fn = str
            filename to save demux artifact
        summary_fn = str
            filename to save demux summary visualization
        qual_plot_fn = str
            filename to save fastq quality plot visualization
        n_qual_plots = int
            number of fastq quality plots to print
    '''

    # import fastq to qiime artifact
    seq_artifact = qiime.Artifact.import_data("RawSequences", seqs_dir)

    # demultiplex
    barcodes = qiime.metadata.MetadataCategory.load(sample_md, index_col)
    demux_seqs = demux.methods.emp(seqs = seq_artifact,
                                   barcodes = barcodes,
                                   rev_comp_barcodes = rc_barcodes,
                                   rev_comp_mapping_barcodes = rc_map_barcodes)
    demux_seqs.per_sample_sequences.save(join(community_dir, demux_fn))

    # demultiplexing summary
    try:
        dsum = demux.visualizers.summarize(demux_seqs.per_sample_sequences)
        dsum.visualization.save(join(community_dir, summary_fn))
    except TypeError:
        # Fails if N=1 https://github.com/qiime2/q2-demux/issues/20
        print("Could not print demux summary: TypeError")

    # view fastq quality plots
    qualplot = dada2.visualizers.plot_qualities(n = n_qual_plots,
                          demultiplexed_seqs = demux_seqs.per_sample_sequences)
    qualplot.visualization.save(join(community_dir, qual_plot_fn))


def denoise_to_phylogeny(communities,
                         mock_data_dir,
                         trim_params,
                         demux_seqs_fn='demux-seqs.qza',
                         rep_seqs_fn='rep_seqs.qza',
                         feature_table_fn='feature_table.qza',
                         summary_fn='feature_table_summary.qzv'):
    '''SampleData[SequencesWithQuality] -> FeatureData[Sequence] +
                                           FeatureTable[Frequency]
        denoise fastqs with dada2, create feature table, rep_seqs,
        and view stats on batch of files.

        communities: LIST
            list of mock communities to extract
        mock_data_dir = filepath
            path to data dir containing communities
        trim_params = dict
            DICT of TUPLES listing dada2 trimming parameters.
            {community : (trim_left, trunc_len)}
        demux_seqs_fn = str
            filename of SampleData[SequencesWithQuality] Artifact to load
        rep_seqs_fn = str
            filename of representative sequences output Artifact
        feature_table_fn = str
            filename of feature table output Artifact
        summary_fn = str
            filename of feature table summary output visualization
    '''
    for community in communities:
        trim_left , trunc_len, buildtree = trim_params[community]
        community_dir = join(mock_data_dir, community)

        # denoise with dada2
        demux_seqs = qiime.Artifact.load(join(community_dir, demux_seqs_fn))
        biom_table, rep_seqs = denoise_to_feature_table(demux_seqs,
                                                        trim_left,
                                                        trunc_len,
                                                        community_dir)
        # Build phylogeny
        if buildtree is True:
            seqs_to_tree(rep_seqs, community_dir)

        print("{0} complete".format(community))


def denoise_to_feature_table(demux_seqs,
                             trim_left,
                             trunc_len,
                             community_dir,
                             rep_seqs_fn='rep_seqs.qza',
                             feature_table_fn='feature_table.qza',
                             summary_fn='feature_table_summary.qzv'):
    '''SampleData[SequencesWithQuality] -> FeatureData[Sequence] +
                                           FeatureTable[Frequency]
    denoise fastqs with dada2, create feature table, rep_seqs,
        and view stats.

        demux_seqs = SampleData[SequencesWithQuality]
            demultiplexed seqs output from qiime.demux.methods.emp()
        trim_left = int
            trim X bases from 5' end
        trunc_len = int
            length to truncate all sequences
        community_dir: path
            destination directory to print results
        rep_seqs_fn = str
            filename of representative sequences output Artifact
        feature_table_fn = str
            filename of feature table output Artifact
        summary_fn = str
            filename of feature table summary output visualization
    '''
    biom_table, rep_seqs = dada2.methods.denoise(demux_seqs,
                                                 trim_left = trim_left,
                                                 trunc_len = trunc_len)
    rep_seqs.save(join(community_dir, rep_seqs_fn))
    biom_table.save(join(community_dir, feature_table_fn))

    # summarize feature table
    feature_table_summary = feature_table.visualizers.summarize(biom_table)
    feature_table_summary.visualization.save(join(community_dir, summary_fn))

    return biom_table, rep_seqs


def seqs_to_tree(rep_seqs, community_dir, filename='phylogeny.qza'):
    '''FeatureData[Sequence] -> phylogeny

    rep_seqs: representative sequences in qiime2 FeatureData[Sequence] Artifact
    community_dir: destination directory to print results
    '''
    # aligned_seqs = alignment.methods.mafft(qiime.Artifact.load(join(\
    #                                      community_dir, 'rep_seqs.qza')))
    aligned_seqs = alignment.methods.mafft(rep_seqs)
    m_aln = alignment.methods.mask(aligned_seqs.alignment)
    unrooted_tree = phylogeny.methods.fasttree(m_aln.masked_alignment)
    tree = phylogeny.methods.midpoint_root(unrooted_tree.tree)
    tree.rooted_tree.save(join(community_dir, filename))
