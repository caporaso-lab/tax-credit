#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2016--, tax-credit development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os import makedirs
from os.path import exists, join, basename
from shutil import copyfile, rmtree
from urllib.request import urlretrieve
from skbio import TreeNode, io
from biom import load_table, Table
from biom.cli.util import write_biom_table
import qiime2
from qiime2.plugins import feature_table, demux, dada2, alignment, phylogeny
from q2_types.feature_data import DNAIterator

from tax_credit.taxa_manipulator import import_taxonomy_to_dict


def extract_mockrobiota_dataset_metadata(mockrobiota_dir, communities):
    '''Extract mock community metadata from mockrobiota dataset-metadata.tsv
    files
    mockrobiota_dir: PATH to mockrobiota directory
    communities: LIST of mock communities to extract
    '''
    dataset_metadata_dict = dict()
    for community in communities:
        dataset_metadata = import_taxonomy_to_dict(
            join(mockrobiota_dir, "data", community, "dataset-metadata.tsv"))
        dataset_metadata_dict[community] = \
            (dataset_metadata['raw-data-url-forward-read'],
             dataset_metadata['raw-data-url-index-read'],
             dataset_metadata['target-gene'])
    return dataset_metadata_dict


def amend_biom_taxonomy_ids(biom_table,
                            empty_taxonomy=['k__', 'p__', 'c__', 'o__',
                                            'f__', 'g__', 's__'],
                            clean_obs_ids=True):
    '''Convert biom table taxonomy strings so that strings with incomplete
    taxonomies are filled out with ambiguous labels
    '''
    if clean_obs_ids is True:
        clean_taxonomy_ids(biom_table)

    new_ids = {}
    for taxa in biom_table.ids(axis='observation'):
        old_taxonomy = taxa.split(';')
        if len(old_taxonomy) < len(empty_taxonomy):
            new_taxonomy = empty_taxonomy
            for i in range(len(old_taxonomy)):
                new_taxonomy[i] = old_taxonomy[i]
            new_ids[taxa] = ';'.join(new_taxonomy)
        else:
            new_ids[taxa] = ';'.join(old_taxonomy)

    return biom_table.update_ids(new_ids, axis='observation')


def clean_taxonomy_ids(table, delete_chars='[]()'):
    new_ids = {obs_id: obs_id.translate(str.maketrans('', '', delete_chars))
               for obs_id in table.ids(axis="observation")}
    return table.update_ids(new_ids, axis='observation')


def extract_mockrobiota_data(communities, community_md, ref_dbs,
                             mockrobiota_dir, mock_data_dir,
                             expected_data_dir, biom_fn='table.L6-taxa.biom'):
    '''Extract sample metadata, raw data files, and expected taxonomy

    from mockrobiota, copy to new destination
    communities: LIST of mock communities to extract
    community_md: DICT of metadata for mock community.
        see extract_mockrobiota_dataset_metadata()
    ref_dbs = DICT mapping marker_gene to reference set names
    mockrobiota_dir = PATH to mockrobiota repo directory
    mock_data_dir = PATH to destination directory
    expected_data_dir = PATH to destination for expected taxonomy files
    '''
    for community in communities:
        # extract dataset metadata/params
        forward_read_url, index_read_url, marker_gene = community_md[community]
        ref_outdir, ref_indir, ref_version, otu_id = ref_dbs[marker_gene][0:4]

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
            if not exists(destination) and file_url_dest[0] != 'NA':
                try:
                    urlretrieve(file_url_dest[0], destination)
                except ValueError:
                    print('Error retrieving {0}'.format(file_url_dest[0]))

        # new directory containing expected taxonomy assignments at each level
        expected_taxa_dir = join(expected_data_dir, community,
                                 ref_outdir, "expected")
        if not exists(expected_taxa_dir):
            makedirs(expected_taxa_dir)
        # copy expected taxonomy.tsv and convert to biom
        exp_taxa_fp = join(expected_taxa_dir, 'expected-taxonomy.tsv')
        exp_biom_fp = join(expected_taxa_dir, biom_fn)
        copyfile(join(mockrobiota_community_dir, ref_indir,
                      ref_version, otu_id, 'expected-taxonomy.tsv'),
                 exp_taxa_fp)
        newbiom = amend_biom_taxonomy_ids(load_table(exp_taxa_fp))
        # add taxonomy ids (names) as observation metadata
        metadata = {sid: {'taxonomy': sid.split(';')}
                    for sid in newbiom.ids(axis='observation')}
        newbiom.add_metadata(metadata, 'observation')
        write_biom_table(newbiom, 'hdf5', exp_biom_fp)


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
        if demux_params[community][0] is True:
            demux_to_plot_quality(seqs_dir=seqs_dir,
                                  sample_md=sample_md,
                                  community_dir=community_dir,
                                  index_col=index_col,
                                  rc_barcodes=demux_params[community][1],
                                  rc_map_barcodes=demux_params[community][2],
                                  demux_fn=demux_fn,
                                  summary_fn=summary_fn,
                                  qual_plot_fn=qual_plot_fn,
                                  n_qual_plots=n_qual_plots)

        else:
            load_demux_seqs(community_dir, seqs_dir, demux_fn, sample_md)

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

    # import fastq to qiime2 artifact
    seq_artifact = qiime2.Artifact.import_data("RawSequences", seqs_dir)

    # demultiplex
    barcodes = qiime2.metadata.MetadataCategory.load(sample_md, index_col)
    demux_seqs = demux.methods.emp_single(
        seqs=seq_artifact, barcodes=barcodes, rev_comp_barcodes=rc_barcodes,
        rev_comp_mapping_barcodes=rc_map_barcodes)
    demux_seqs.per_sample_sequences.save(join(community_dir, demux_fn))

    visualize_qual(demux_seqs, community_dir, summary_fn,
                   qual_plot_fn, n_qual_plots)


def visualize_qual(demux_seqs, community_dir, summary_fn,
                   qual_plot_fn, n_qual_plots=1):

    '''visualize demux summary and fastq quality plots.'''

    # demultiplexing summary
    try:
        dsum = demux.visualizers.summarize(demux_seqs.per_sample_sequences)
        dsum.visualization.save(join(community_dir, summary_fn))
    except TypeError:
        # Fails if N=1 https://github.com/qiime2/q2-demux/issues/20
        print("Could not print demux summary: TypeError")

    # view fastq quality plots
    qualplot = dada2.visualizers.plot_qualities(
        n=n_qual_plots, demultiplexed_seqs=demux_seqs.per_sample_sequences)
    qualplot.visualization.save(join(community_dir, qual_plot_fn))


def load_demux_seqs(community_dir, seqs_dir, demux_fn, sample_md):
    '''Load demultiplexed sequences into artifact and write artifact.'''

    # extract sample name
    with open(sample_md, 'r') as md_in:
        md_in.readline()
        sample_id = md_in.readline().strip().split('\t')[0]

    # assemble artifact components
    tmpdir = join(seqs_dir, 'tmp')
    seq_fn = '{0}_1_L001_R1_001.fastq.gz'.format(sample_id)
    if not exists(tmpdir):
        makedirs(tmpdir)
    copyfile(join(seqs_dir, 'sequences.fastq.gz'),
             join(tmpdir, seq_fn))
    with open(join(tmpdir, 'MANIFEST'), 'w') as manifest:
        manifest.write('sample-id,filename,direction\n')
        manifest.write('{0},{1},forward'.format(sample_id, seq_fn))
    with open(join(tmpdir, 'metadata.yml'), 'w') as yml:
        yml.write('\n')

    seqs = qiime2.Artifact.import_data(
        "SampleData[SequencesWithQuality]", tmpdir)
    rmtree(tmpdir)

    seqs.save(join(community_dir, demux_fn))

    # Visualize currently fails on artifacts containing single sample
    # (seqs has no attribute 'per_sample_sequences')
    # visualize_qual(seqs, community_dir, summary_fn,
    #               qual_plot_fn, n_qual_plots)


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
        trim_left, trunc_len, buildtree = trim_params[community]
        community_dir = join(mock_data_dir, community)

        # denoise with dada2
        demux_seqs = qiime2.Artifact.load(join(community_dir, demux_seqs_fn))
        biom_table, rep_seqs = denoise_to_feature_table(
            demux_seqs, trim_left, trunc_len, community_dir)

        # Build phylogeny
        if buildtree is True:
            seqs_to_tree(rep_seqs, community_dir)

        print("{0} complete".format(community))


def denoise_to_feature_table(demux_seqs,
                             trim_left,
                             trunc_len,
                             community_dir,
                             rep_seqs_fn='rep_seqs',
                             feature_table_fn='feature_table.qza',
                             biom_table_fn='feature_table.biom',
                             summary_fn='feature_table_summary.qzv'):
    '''SampleData[SequencesWithQuality] -> FeatureData[Sequence] +
                                           FeatureTable[Frequency]
    denoise fastqs with dada2, create feature table, rep_seqs,
        and view stats.

        demux_seqs = SampleData[SequencesWithQuality]
            demultiplexed seqs output from qiime2.demux.methods.emp()
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
    biom_table, rep_seqs = dada2.methods.denoise_single(
        demux_seqs, trim_left=trim_left, trunc_len=trunc_len)
    # save Artifact
    rep_seqs.save(join(community_dir, rep_seqs_fn))

    # save biom Artifact
    biom_table.save(join(community_dir, feature_table_fn))
    biom_table_fp = join(community_dir, biom_table_fn)
    write_biom_table(biom_table.view(Table), 'hdf5', biom_table_fp)

    # summarize feature table
    feature_table_summary = feature_table.visualizers.summarize(biom_table)
    feature_table_summary.visualization.save(join(community_dir, summary_fn))

    return biom_table, rep_seqs


def seqs_to_tree(rep_seqs, community_dir, filename='phylogeny.qza'):
    '''FeatureData[Sequence] -> phylogeny

    rep_seqs: FeatureData[Sequence] Artifact
        representative sequences from dada2
    community_dir: path
        destination directory to print results
    '''
    aligned_seqs = alignment.methods.mafft(rep_seqs)
    m_aln = alignment.methods.mask(aligned_seqs.alignment)
    unrooted_tree = phylogeny.methods.fasttree(m_aln.masked_alignment)
    tree = phylogeny.methods.midpoint_root(unrooted_tree.tree)
    tree.rooted_tree.save(join(community_dir, filename))


def transport_to_repo(communities,
                      mock_data_dir,
                      project_dir,
                      sample_type_dirname='mock-community',
                      rep_seqs_fn='rep_seqs.qza',
                      feature_table_fn='feature_table.qza',
                      tree_fn='phylogeny.qza',
                      sample_md_fn='sample-metadata.tsv',
                      biom_table_fn='feature_table.biom',
                      fasta_fn='rep_seqs.fna',
                      newick_fn='phylogeny.tre'):
    '''Copy essential mock community data to tax-credit repo

    communities: list
        list of dir names in mock_data_dir, a.k.a. names of mock communities
    mock_data_dir: path
        source directory containing mock communities dirs of results
    project_dir: path
        path to tax-credit repo directory
    sample_type_dirname: str
        name of destination directory to contain communities dirs. The analog
        of mock_data_dir in the repo, dirs for individual communities will be
        located in project_dir/data/sample_type_dirname/community
    rep_seqs_fn: str
        name of rep seqs FeatureData[Sequence] Artifact in community_dir
    feature_table_fn: str
        name of rep seqs FeatureTable[Frequency] Artifact in community_dir
    tree_fn: str
        name of Phylogeny[Rooted] Artifact in community_dir
    sample_md_fn: str
        name of metadata mapping file in community_dir
    biom_table_fn: str
        destination name of biom table in project_dir
    fasta_fn: str
        destination name of fasta file in project_dir
    newick_fn: str
        destination name of newick format tree in project_dir
    '''

    for community in communities:
        community_dir = join(mock_data_dir, community)

        # Define base dir destination for mock community directories
        repo_destination = join(project_dir, "data",
                                sample_type_dirname, community)
        if not exists(repo_destination):
            makedirs(repo_destination)

        # Files to move
        rep_seqs = join(community_dir, rep_seqs_fn)
        feature_table = join(community_dir, feature_table_fn)
        tree = join(community_dir, tree_fn)
        sample_md = join(community_dir, sample_md_fn)

        biom_table_fp = join(community_dir, biom_table_fn)
        rep_seqs_fp = join(community_dir, fasta_fn)
        tree_fp = join(community_dir, newick_fn)

        # Extract biom, tree, rep_seqs
        rep_seqs_fna = qiime2.Artifact.load(rep_seqs).view(DNAIterator)
        io.write(rep_seqs_fna.generator, format='fasta', into=rep_seqs_fp)

        if exists(tree):
            qiime2.Artifact.load(tree).view(TreeNode).write(tree_fp)

        # Move to repo:
        for f in [rep_seqs, feature_table, tree, sample_md,
                  biom_table_fp, rep_seqs_fp, tree_fp]:
            if exists(f):
                copyfile(f, join(repo_destination, basename(f)))
