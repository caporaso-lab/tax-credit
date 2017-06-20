#!/usr/bin/env python


# ----------------------------------------------------------------------------
# Copyright (c) 2016--, tax-credit development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


import os
import glob
import csv
import shutil
import re
import skbio
import qiime2
import pandas
import seaborn as sns
import matplotlib
import numpy


def ref_db_to_dict(ref_dbs):
    '''Turn reference database into a dictionary.

    ref_dbs: list of tuples
        [(name, fasta path, taxonomy path)]

    '''
    refs = {}
    taxs = {}
    for ref_name, db_file, tax_file in ref_dbs:
        seqs = skbio.io.read(db_file, format='fasta', constructor=skbio.DNA)
        refs[ref_name] = list(seqs)
        tax_map = {}
        with open(tax_file) as tax_fh:
            reader = csv.reader(tax_fh, delimiter='\t')
            for row in reader:
                tax_map[row[0]] = row[1]
            taxs[ref_name] = tax_map
    return refs, taxs


def single_character_matches(char, query, subject):
    query = ((query == char) + (query == 'N')).astype(int)
    subject = ((subject == char) + (subject == 'N')).astype(int)
    return numpy.convolve(query[::-1], subject)


def gapless_mismatches(query, subject):
    matches = single_character_matches('A', query, subject)
    matches += single_character_matches('C', query, subject)
    matches += single_character_matches('G', query, subject)
    matches += single_character_matches('T', query, subject)
    return len(query) - matches.max()


def compare_observed_and_expected(observed, expected):
    ids = [obs.metadata['id'] for obs in observed]
    rev_observed = [numpy.array(list(str(obs.reverse_complement())))
                    for obs in observed]
    observed = [numpy.array(list(str(obs))) for obs in observed]
    expected = list(expected)
    keys = [' '.join([exp.metadata['id'], exp.metadata['description']]).strip()
            for exp in expected]
    expected = [numpy.array(list(str(exp))) for exp in expected]
    comparison = {}
    for _id, obs, rev in zip(ids, observed, rev_observed):
        matches = [(gapless_mismatches(obs, exp), key) for
                   exp, key in zip(expected, keys)]
        matches += [(gapless_mismatches(rev, exp), key) for
                    exp, key in zip(expected, keys)]
        matches.sort()
        comparison[_id] = [(e, k) for e, k in matches if e == matches[0][0]]
    return comparison


def match_expected_seqs_to_taxonomy(data_dir, mockrobiota_dir, refs, taxs):
    '''Map expected_seqs to taxonomy assignment matched to reference database
    of interest.

    data_dir: path
        tax-credit repository base directory.
    mockrobiota_dir: path
        mockrobiota repository base directory.
    refs: dict of list
        {ref_name: [ref sequences]}
    taxs: dict of dicts
        {ref_name: {id_: taxonomy}}
    '''

    # *** This function is way too long! We should break it up!

    # *** We also need to make the "save the taxonomies of the matched reads"
    # *** better and more automated, e.g., to find all ref dbs in
    # *** com_dir and match to ref taxonomies in ref_dbs. This will also
    # *** require searching for appropriate taxonomies in this_mockrobiota_dir

    data_dir = os.path.join(data_dir, 'data')
    mockrobiota_dir = os.path.join(mockrobiota_dir, 'data')
    mock_dir = os.path.join(data_dir, 'mock-community')
    expected_dir = os.path.join(
        data_dir, 'precomputed-results', 'mock-community')

    for com_dir in glob.glob(os.path.join(expected_dir, '*')):
        if not os.path.isdir(com_dir):
            continue
        community = os.path.basename(com_dir)
        this_mock_dir = os.path.join(mock_dir, community)
        this_mockrobiota_dir = os.path.join(mockrobiota_dir, community)

        # find the expected taxonomies
        if os.path.exists(os.path.join(this_mockrobiota_dir, 'greengenes')):
            identifier_file_dir = os.path.join(
                this_mockrobiota_dir, 'greengenes', '13-8', '99-otus')
            ref_seqs = refs['greengenes']
            taxonomy_map = taxs['greengenes']
        elif os.path.exists(os.path.join(this_mockrobiota_dir, 'unite')):
            identifier_file_dir = os.path.join(
                this_mockrobiota_dir, 'unite', '7-1', '99-otus')
            ref_seqs = refs['unite']
            taxonomy_map = taxs['unite']
        else:
            raise RuntimeError('could not find identifiers for ' + community)

        # set the destinations
        es_fp = os.path.join(this_mock_dir, 'expected-sequences.fasta')
        et_fp = os.path.join(this_mock_dir, 'matched-sequence-taxonomies.tsv')
        identifier_file = os.path.join(
            identifier_file_dir, 'database-identifiers.tsv')

        # check to see whether the expected sequences are precompiled ...
        expected_expected = os.path.join(
            this_mockrobiota_dir, 'source', 'expected-sequences.fasta')
        if os.path.exists(expected_expected):
            print(community)
            expected_seqs = skbio.io.read(
                expected_expected, format='fasta', constructor=skbio.DNA)
            expected_ids = {' '.join([s.metadata['id'],
                            s.metadata['description']]).strip()
                            for s in expected_seqs}
            shutil.copy(expected_expected, es_fp)
            est_fp = os.path.join(
                this_mock_dir, 'expected-sequence-taxonomies.tsv')
            with open(est_fp) as est_fh:
                with open(et_fp, 'w') as et_fh:
                    reader = csv.reader(est_fh, delimiter='\t')
                    writer = csv.writer(et_fh, delimiter='\t')
                    writer.writerow(next(reader))
                    for tax_id, tax in reader:
                        est_id = tax_id.split(';')[-1]
                        est_ids = re.split('_| ', est_id)
                        # attempt to find exact match
                        if est_id in expected_ids:
                            writer.writerow([est_id, tax])
                        # else find near match in expected_ids
                        else:
                            for exp_id in expected_ids:
                                # instead of hard-coding rules for individual
                                # datasets let's just do a simple search with
                                # flexible delimiters attempt match in full
                                # name
                                if (re.search(est_id, exp_id) or
                                        # last resort: just match genus/species
                                        (re.search(est_ids[0], exp_id) and
                                         re.search(est_ids[1], exp_id))):
                                    writer.writerow([exp_id, tax])
                                    # break

        # ... and otherwise derive the expected sequences from the identifiers
        # and references
        elif os.path.exists(identifier_file):
            with open(identifier_file) as csvfile:
                identifiers = csv.reader(csvfile, delimiter='\t')
                expected_seqs = []
                for row in identifiers:
                    for _id in row[1].split():
                        for seq in ref_seqs:
                            if _id == seq.metadata['id']:
                                expected_seqs.append(seq)
                                break
                        else:
                            raise RuntimeError('could not find ' + _id +
                                               ' for ' + community)

            # now write them out for safe-keeping
            skbio.io.write((s for s in expected_seqs), 'fasta', es_fp)
            print(community)
            with open(et_fp, 'w') as et_fh:
                writer = csv.writer(et_fh, delimiter='\t')
                writer.writerow(['Taxonomy', 'Standard Taxonomy'])
                for seq in expected_seqs:
                    _id = seq.metadata['id']
                    writer.writerow([_id, taxonomy_map[_id]])
        else:
            print('skipping ' + community)
            print(expected_expected)


def generate_trueish_taxonomies(data_dir, force=False):
    '''Generate "true" taxonomy assignment for each observed sequence in mock
    community, and plot mismatch distributions across all observed sequences
    in each mock community.

    data_dir: path
        tax-credit repository base directory.
    '''
    # *** We should split this up into at least two functions: generate
    # *** trueish taxonomies and plot mismatch distributions.

    data_dir = os.path.join(data_dir, 'data')
    mock_dir = os.path.join(data_dir, 'mock-community')
    expected_dir = os.path.join(
        data_dir, 'precomputed-results', 'mock-community')

    for com_dir in glob.glob(os.path.join(expected_dir, '*')):
        if not force and not os.path.isdir(com_dir):
            continue
        community = os.path.basename(com_dir)
        this_mock_dir = os.path.join(mock_dir, community)

        # load the taxonomy map
        tax_file = os.path.join(
            this_mock_dir, 'matched-sequence-taxonomies.tsv')
        if os.path.exists(tax_file):
            with open(tax_file) as tf_fh:
                reader = csv.reader(tf_fh, delimiter='\t')
                next(reader)
                tax_map = {seq: tax for seq, tax in reader}

            # load the expected sequences
            expected_file = os.path.join(
                this_mock_dir, 'expected-sequences.fasta')
            expected_sequences = skbio.io.read(
                expected_file, format='fasta', constructor=skbio.DNA)

            # compare with the observed sequences
            observed_file = os.path.join(this_mock_dir, 'rep_seqs.fna')
            observed = skbio.io.read(
                observed_file, format='fasta', constructor=skbio.DNA)
            observed = list(observed)
            mismatches = compare_observed_and_expected(
                observed, expected_sequences)

            # save the taxonomies of the matched reads
            # we need them in the "expected results" directories
            print(community)
            if os.path.exists(os.path.join(
                    com_dir, 'gg_13_8_otus', 'expected')):
                result_file = os.path.join(com_dir, 'gg_13_8_otus', 'expected',
                                           'trueish-taxonomies.tsv')
            elif os.path.exists(os.path.join(
                    com_dir, 'unite_20.11.2016_clean_fullITS', 'expected')):
                result_file = os.path.join(
                    com_dir, 'unite_20.11.2016_clean_fullITS', 'expected',
                    'trueish-taxonomies.tsv')
            else:
                raise RuntimeError('Must specify expected results path for '
                                   '{0}'.format(com_dir))

            max_mismatches = 3
            with open(result_file, 'w') as result_fh:
                writer = csv.writer(result_fh, delimiter='\t')
                writer.writerow(['id', 'taxonomy'])
                for rep_id, mismatches_for_id in mismatches.items():
                    if mismatches_for_id[0][0] > max_mismatches:
                        writer.writerow([rep_id, 'other'])
                        continue
                    matched_taxonomies = set()
                    for m, e in mismatches_for_id:
                        if e in tax_map:
                            matched_taxonomies.add(tax_map[e])
                        else:
                            matched_taxonomies.add('NO-EXPECTED-TAXONOMY')
                            print('WARNING: no expected taxonomy for ' + e)
                    matched_taxonomies = list(matched_taxonomies)
                    if len(matched_taxonomies) > 1:
                        print('WARNING: ' + rep_id + ' matches')
                        print(' and\n'.join(matched_taxonomies))
                        print('with ' + str(mismatches_for_id[0][0]) +
                              ' mismatches')
                    for matched_taxonomy in matched_taxonomies:
                        writer.writerow([rep_id, matched_taxonomy])

            feature_file = os.path.join(this_mock_dir, 'feature_table.qza')
            features = qiime2.Artifact.load(feature_file)
            features = features.view(pandas.DataFrame)
            best_mismatches = list(mm[0][0] for mm in mismatches.values())
            weighted_mismatches = []
            found = 0
            total = 0
            for obs in mismatches:
                try:
                    count = int(features[obs])
                except TypeError:
                    count = int(sum(features[obs]))
                num_mismatches = mismatches[obs][0][0]
                weighted_mismatches.extend([num_mismatches]*count)
                if num_mismatches <= max_mismatches:
                    found += count
                total += count
            print('Guessed taxonomy for %d of %d reads (%.1f%%)' %
                  (found, total, 100*found/total))
            sns.distplot(best_mismatches, kde=False, bins=50,
                         axlabel='Number of Mismatches')
            matplotlib.pyplot.show()
            sns.distplot(weighted_mismatches, kde=False, bins=50,
                         axlabel='Number of Mismatches')
            matplotlib.pyplot.show()
