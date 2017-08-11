#!/usr/bin/env python


# ----------------------------------------------------------------------------
# Copyright (c) 2016--, tax-credit development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import join, exists, split, sep, expandvars, basename, splitext
from os import makedirs, remove, system, stat
from glob import glob
from itertools import product
from shutil import rmtree, move
from random import sample
from time import time
from collections import Counter, defaultdict

import matplotlib.pyplot as plt
from seaborn import regplot
from biom import load_table
from biom.cli.util import write_biom_table
from biom.parse import MetadataMap
from numpy import random
from skbio import io
import pandas as pd
from sklearn.model_selection import StratifiedKFold

from tax_credit.taxa_manipulator import (accept_list_or_file,
                                         import_to_list,
                                         import_taxonomy_to_dict,
                                         export_list_to_file,
                                         filter_sequences,
                                         string_search,
                                         trim_taxonomy_strings,
                                         extract_taxa_names)
from qiime2 import Artifact
from qiime2.plugins import feature_classifier
from q2_types.feature_data import DNAIterator


def gen_param_sweep(data_dir, results_dir, reference_dbs,
                    dataset_reference_combinations,
                    method_parameters_combinations,
                    output_name='rep_set_tax_assignments.txt', force=False):
    '''Create list of commands from input dictionaries of method_parameter and
    dataset_reference_combinations
    '''
    for dataset, reference in dataset_reference_combinations:
        input_dir = join(data_dir, dataset)
        reference_seqs, reference_tax = reference_dbs[reference]
        for method, parameters in method_parameters_combinations.items():
            parameter_ids = sorted(parameters.keys())
            for parameter_combination in product(*[parameters[id_]
                                                 for id_ in parameter_ids]):
                parameter_comb_id = ':'.join(map(str, parameter_combination))
                parameter_output_dir = join(results_dir, dataset, reference,
                                            method, parameter_comb_id)
                if not exists(join(parameter_output_dir, output_name)) \
                        or force:
                    params = dict(zip(parameter_ids, parameter_combination))
                    yield (parameter_output_dir, input_dir, reference_seqs,
                           reference_tax, method, params)


def parameter_sweep(data_dir, results_dir, reference_dbs,
                    dataset_reference_combinations,
                    method_parameters_combinations,
                    command_template="mkdir -p {0} ; assign_taxonomy.py -v -i "
                                     "{1} -o {0} -r {2} -t {3} -m {4} {5} "
                                     "--rdp_max_memory 16000",
                    infile='rep_set.fna',
                    output_name='rep_set_tax_assignments.txt', force=False):
    '''Create list of commands from input dictionaries of method_parameter and
    dataset_reference_combinations
    '''
    sweep = gen_param_sweep(data_dir, results_dir, reference_dbs,
                            dataset_reference_combinations,
                            method_parameters_combinations,
                            output_name, force)
    commands = []
    for assignment in sweep:
        (parameter_output_dir, input_dir, reference_seqs, reference_tax,
         method, params) = assignment
        query_seqs = join(input_dir, infile)
        parameter_str = ' '.join(['--%s %s' % e for e in params.items()])
        command = command_template.format(parameter_output_dir, query_seqs,
                                          reference_seqs, reference_tax,
                                          method, parameter_str)
        commands.append(command)
    return commands


def add_metadata_to_biom_table(biom_input_fp, taxonomy_map_fp, biom_output_fp):
    '''Load biom, add metadata, write to new table'''
    newbiom = load_table(biom_input_fp)
    if stat(taxonomy_map_fp).st_size == 0:
        metadata = {}
    else:
        metadata = MetadataMap.from_file(
            taxonomy_map_fp, header=['Sample ID', 'taxonomy', 'c'])
    newbiom.add_metadata(metadata, 'observation')
    write_biom_table(newbiom, 'json', biom_output_fp)


def generate_per_method_biom_tables(taxonomy_glob, data_dir,
                                    biom_input_fn='feature_table.biom',
                                    biom_output_fn='table.biom'):
    '''Create biom tables from glob of taxonomy assignment results'''
    taxonomy_map_fps = glob(expandvars(taxonomy_glob))
    for taxonomy_map_fp in taxonomy_map_fps:
        dataset_id = taxonomy_map_fp.split(sep)[-5]
        biom_input_fp = join(data_dir, dataset_id, biom_input_fn)
        output_dir = split(taxonomy_map_fp)[0]
        biom_output_fp = join(output_dir, biom_output_fn)
        if exists(biom_output_fp):
            remove(biom_output_fp)
        add_metadata_to_biom_table(biom_input_fp, taxonomy_map_fp,
                                   biom_output_fp)


def move_results_to_repository(method_dirs, precomputed_results_dir):
    '''Move new taxonomy assignment results to git repository'''
    for method_dir in method_dirs:
        f = method_dir.split(sep)
        dataset_id, db_id, method_id, param_id = f[-4], f[-3], f[-2], f[-1]

        new_location = join(precomputed_results_dir, dataset_id, db_id,
                            method_id, param_id)
        if exists(new_location):
            rmtree(new_location)
        move(method_dir, new_location)


def clean_database(taxa_in, seqs_in, db_dir,
                   junk='__;|__$|_sp$|unknown|unidentified'):

    '''Remove ambiguous and empty taxonomies from reference seqs/taxonomy.

    taxa_in: path
        File containing taxonomy strings in tab-separated format:
        <SequenceID>    <taxonomy string>

    seqs_in: path
        File containing sequences corresponding to taxa_in, in fasta format.

    db_dir: dir path
        Output directory.

    junk: str
        '|'-separated list of search terms. Taxonomies containing these terms
        will be removed from the database.
    '''

    clean_taxa = join(
        db_dir, '{0}_clean.tsv'.format(basename(splitext(taxa_in)[0])))
    clean_fasta = join(
        db_dir, '{0}_clean.fasta'.format(basename(splitext(seqs_in)[0])))

    # Remove empty taxa from ref taxonomy
    taxa = string_search(taxa_in, junk, discard=True)
    # Remove brackets (and other special characters causing problems)
    clean_list = [line.translate(str.maketrans('', '', '[]()'))
                  for line in taxa]
    export_list_to_file(clean_list, clean_taxa)
    # Remove corresponding seqs from ref fasta
    filter_sequences(seqs_in, clean_fasta, clean_taxa, keep=True)

    return clean_taxa, clean_fasta


def seq_count(infile):
    '''Count sequences in fasta file'''
    with open(infile, "r") as f:
        count = sum(1 for line in f) / 2
    return count


class Node(object):
    def __init__(self):
        self.tips = set()
        self.children = defaultdict(Node)

    def add(self, sid, taxon):
        self.tips.add(sid)
        if taxon:
            self.children[taxon[0]].add(sid, taxon[1:])

    def __lt__(self, other):
        return len(self.tips) < len(other.tips)


def build_tree(taxonomy, sequences):
    tree = Node()
    filtered_ids = {s.metadata['id'] for s in sequences}
    for sid, taxon in taxonomy.view(pd.Series).items():
        if sid not in filtered_ids:
            continue
        tree.add(sid, taxon.split(';'))
    return tree


def get_strata(tree, k, taxon=[]):
    if not tree.children:
        return [(';'.join(taxon), tree.tips)]
    sorted_children = map(list, map(reversed, tree.children.items()))
    sorted_children = iter(sorted(sorted_children))
    misc = set()
    # aggregate children with small tip sets
    for child, level in sorted_children:
        if len(child.tips) >= k:
            break
        misc.update(child.tips)
    else:  # all the tips are in misc
        return [(';'.join(taxon), misc)]
    # grab the tips from the child that's not in misc
    strata = get_strata(child, k, taxon+[level])
    # get the rest
    for child, level in sorted_children:
        strata.extend(get_strata(child, k, taxon+[level]))
    # if there were more than k misc, make a stratum for them
    if len(misc) > k:
        strata = [(';'.join(taxon), misc)] + strata
    else:  # randomly add them to the other strata
        for i, tip in enumerate(misc):
            strata[i % len(strata)][1].add(tip)
    return strata


def distances(A, B, n):
    iA = random.choice(range(len(A)), size=n)
    iB = random.choice(range(len(B)), size=n)
    d = [0]*n
    for i, (a, b) in enumerate(zip(iA, iB)):
        d[i] = sum(x != y for x, y in zip(A[a], B[b]))
    return d


def distance_comparison(dataframe, data_dir, test_name, samples=10000):
    simulated_dir = join(data_dir, test_name)
    for index, data in dataframe.iterrows():
        lengths = Counter()
        inner = []
        outer = []
        trainsets = glob(join(simulated_dir, index+'*', 'ref_seqs.fasta'))
        testsets = glob(join(simulated_dir, index+'*', 'query.fasta'))
        for train_fp, test_fp in zip(trainsets, testsets):
            train = list(map(str, io.read(train_fp, format='fasta')))
            test = list(map(str, io.read(test_fp, format='fasta')))
            if not train or not test:
                continue
            lengths.update(map(len, test))
            inner.extend(distances(train, train, samples))
            outer.extend(distances(train, test, samples))
        inner.sort()
        outer.sort()
        df = pd.DataFrame(
            {'train/train': inner, 'train/test': outer})

        plt.figure()  # figsize=(width, height))
        ax = regplot('train/train', 'train/test', df, fit_reg=False)
        ax.set_title(index, fontsize=20)
        maxval = max((inner[-1], outer[-1]))
        plt.plot([0, maxval], [0, maxval], linewidth=2)
        plt.show()


def generate_simulated_datasets(dataframe, data_dir, read_length, iterations,
                                levelrange=range(6, 0, -1), force=False,
                                min_read_length=80):
    '''From a dataframe of sequence reference databases, build training/test
    sets of "novel taxa" queries, taxonomies, and reference seqs/taxonomies.

    dataframe: pandas.DataFrame
    '''
    if iterations < 2:
        raise ValueError('Must perform two or more iterations for '
                         'construction of cross-validated datasets.')

    cv_dir = join(data_dir, 'cross-validated')
    novel_dir = join(data_dir, 'novel-taxa-simulations')
    for index, data in dataframe.iterrows():
        db_dir = join(data_dir, 'ref_dbs', data['Reference id'])
        if not exists(db_dir):
            makedirs(db_dir)

        # Clean taxonomy/sequences, to remove empty/ambiguous taxonomies
        clean_fasta = join(db_dir, '{0}_clean.fasta'.format(
            basename(splitext(data['Reference file path'])[0])))
        clean_taxa = join(db_dir, '{0}_clean.tsv'.format(
            basename(splitext(data['Reference tax path'])[0])))
        if not exists(clean_fasta) or force:
            clean_taxa, clean_fasta = clean_database(
              data['Reference tax path'], data['Reference file path'], db_dir)

        # Extract amplicons and filter by minimum length
        primer_pair = '{0}-{1}'.format(data['Fwd primer id'],
                                       data['Rev primer id'])
        base, ext = splitext(clean_fasta)
        simulated_reads_fp = join('{0}_{1}_trim{2}{3}'.format(
            base, primer_pair, read_length, ext))
        if not exists(simulated_reads_fp) or force:
            # Trim reference sequences to amplicon target
            seqs = Artifact.import_data('FeatureData[Sequence]', clean_fasta)
            trimmed = feature_classifier.methods.extract_reads(
                sequences=seqs, trunc_len=read_length,
                f_primer=data['Fwd primer'], r_primer=data['Rev primer']).reads

            with open(simulated_reads_fp, 'w') as simulated_reads:
                for seq in trimmed.view(DNAIterator):
                    if len(seq) >= min_read_length:
                        seq.write(simulated_reads, format='fasta')
        else:
            print('simulated reads and amplicons exist: skipping extraction')

        # Print sequence counts to see how many seqs were filtered
        print(index, 'Sequence Counts')
        print('Raw Fasta:           ', seq_count(data['Reference file path']))
        print('Clean Fasta:         ', seq_count(clean_fasta))
        print('Simulated Reads:     ', seq_count(simulated_reads_fp))

        # Generate simulated community query and reference seqs/taxa pairs
        generate_cross_validated_sequences(
            clean_taxa, simulated_reads_fp, index, iterations, cv_dir)

        # Generate novel query and reference seqs/taxa pairs
        generate_novel_sequence_sets(cv_dir, novel_dir, levelrange=levelrange)


def generate_novel_sequence_sets(cv_dir, novel_dir,
                                 levelrange=range(6, 0, -1)):
    '''Generate paired query/reference fastas and taxonomies for novel taxa
    analysis, given an input of simulated amplicon taxonomies (read_taxa)
    and fastas (simulated_reads_fp), the index (database) name, # of
    iterations to perform (cross-validated data subsets), the output dir
    (novel_dir).
    novel_dir: path
        base output directory to contain simulated datasets
    cv_dir: path
        directory containing cross-validation data sets for filtering
    '''

    for cv_fold_dir in glob(join(cv_dir, '*')):
        index, REF, iteration = basename(cv_fold_dir).split('-')
        for level in levelrange:
            novel_fold_dir = join(
                novel_dir, '-'.join([index, REF, 'L'+str(level), iteration]))
            if not exists(novel_fold_dir):
                makedirs(novel_fold_dir)

            # Trim taxonomy strings to level X and copy them into the directory
            query_taxa = trim_taxonomy_strings(
                join(cv_fold_dir, 'query_taxa.tsv'), level)
            query_taxa_fp = join(novel_fold_dir, 'query_taxa.tsv')
            export_list_to_file(query_taxa, query_taxa_fp)

            # Create REF TAXONOMY from list of taxonomies that
            #    DO NOT match QUERY taxonomies
            ref_taxa_fp = join(novel_fold_dir, 'ref_taxa.tsv')
            name_list = '^' + '$|^'.join(list(set(extract_taxa_names(
                query_taxa_fp, slice(1, level+1), stripchars='')))) + '$'
            ref_taxa = string_search(join(cv_fold_dir, 'ref_taxa.tsv'),
                                     name_list, discard=True,
                                     field=slice(1, level+1))
            export_list_to_file(ref_taxa, ref_taxa_fp)

            # Now trim the taxonomy strings to level X-1 for consistency with
            # CV and remove any that are unclassified at that level
            query_taxa = trim_taxonomy_strings(query_taxa_fp, level-1)
            separated = defaultdict(list)
            for s in query_taxa:
                separated[s.count(';')].append(s)
            export_list_to_file(separated[max(separated)], query_taxa_fp)

            # Create REF: Filter ref database to contain only seqs that
            #    match non-matching taxonomy strings
            ref_fp = join(novel_fold_dir, 'ref_seqs.fasta')
            filter_sequences(join(cv_fold_dir, 'ref_seqs.fasta'),
                             ref_fp, ref_taxa_fp, keep=True)

            # Create QUERY: Filter ref database to contain only seqs
            #    that match QUERY TAXONOMY
            query_fp = join(novel_fold_dir, 'query.fasta')
            filter_sequences(join(cv_fold_dir, 'query.fasta'),
                             query_fp, query_taxa_fp, keep=True)

            # Encode as Artifacts for convenience
            artifact = Artifact.import_data('FeatureData[Sequence]', ref_fp)
            artifact.save(ref_fp[:-5] + 'qza')
            artifact = Artifact.import_data('FeatureData[Sequence]', query_fp)
            artifact.save(query_fp[:-5] + 'qza')
            artifact = Artifact.import_data(
                'FeatureData[Taxonomy]', ref_taxa_fp,
                view_type='HeaderlessTSVTaxonomyFormat')
            artifact.save(ref_taxa_fp[:-3] + 'qza')


def generate_cross_validated_sequences(read_taxa, simulated_reads_fp, index,
                                       iterations, cv_dir):
    '''Generates simulated community files (fasta and taxonomy) as subsets of
    simulated amplicons/taxa for cross-validated taxonomy assignment. Selects
    duplicated taxa names, evenly allocates these among subsets as query taxa
    (test set), generates ref taxa (training set) that do not match query fasta
    IDs, and creates fasta files to match each of these sets.
    read_taxa: list or path
        list or file of taxonomies corresponding to simulated_reads_fp
    simulated_reads_fp: path
        simulated amplicon reads (fasta format file)
    index: str
        reference database name
    iterations: int >= 2
        number of subsets to create
    cv_dir: path
        base output directory to contain simulated datasets
    '''
    if iterations < 2:
        raise ValueError('Must perform two or more iterations for '
                         'construction of cross-validated datasets.')

    # Stratify the data and form the CV data sets
    simulated_reads = list(io.read(simulated_reads_fp, format='fasta'))
    taxonomy = Artifact.import_data(
        'FeatureData[Taxonomy]', read_taxa,
        view_type='HeaderlessTSVTaxonomyFormat')
    tree = build_tree(taxonomy, simulated_reads)
    strata = get_strata(tree, iterations)
    print(index + ': generating', iterations, 'folds on', len(strata),
          'strata')
    X, y = zip(*[(s, t) for t, ss in strata for s in ss])
    skf = StratifiedKFold(
        n_splits=iterations, shuffle=True, random_state=0)
    splits = []
    for train, test in skf.split(X, y):
        train_set = {X[i] for i in train}
        test_set = {X[i] for i in test}
        splits.append((train_set, test_set))

    # Output the CV data sets in the expected formats
    taxonomy_series = taxonomy.view(pd.Series)
    for iteration, (train, test) in enumerate(splits):
        db_iter_dir = join(cv_dir, '{0}-iter{1}'.format(index, iteration))
        if not exists(db_iter_dir):
            makedirs(db_iter_dir)
        query_taxa_fp = join(db_iter_dir, 'query_taxa.tsv')
        query_fp = join(db_iter_dir, 'query.fasta')
        ref_fp = join(db_iter_dir, 'ref_seqs.fasta')
        ref_taxa_fp = join(db_iter_dir, 'ref_taxa.tsv')

        # Output the taxa files
        train_series = taxonomy_series[train]
        train_series.to_csv(ref_taxa_fp, sep='\t')
        # If a taxonomy in the test set doesn't exist in the training set, trim
        # it until it does
        train_taxonomies = set()
        for taxonomy in train_series.values:
            taxonomy = taxonomy.split(';')
            for level in range(1, len(taxonomy)+1):
                train_taxonomies.add(';'.join(taxonomy[:level]))
        test_list = []
        for sid in test:
            taxonomy = taxonomy_series[sid].split(';')
            for level in range(len(taxonomy), 0, -1):
                if ';'.join(taxonomy[:level]) in train_taxonomies:
                    test_list.append(
                        '\t'.join([sid, ';'.join(taxonomy[:level]).strip()]))
                    break
            else:
                raise RuntimeError('unknown kingdom in query set')
        export_list_to_file(test_list, query_taxa_fp)
        # Output the reference files
        with open(ref_fp, 'w') as ref_fasta:
            with open(query_fp, 'w') as query_fasta:
                for seq in simulated_reads:
                    if seq.metadata['id'] in train:
                        seq.write(ref_fasta, format='fasta')
                    else:
                        seq.write(query_fasta, format='fasta')

        # Encode as Artifacts for convenience
        artifact = Artifact.import_data('FeatureData[Sequence]', ref_fp)
        artifact.save(ref_fp[:-5] + 'qza')
        artifact = Artifact.import_data('FeatureData[Sequence]', query_fp)
        artifact.save(query_fp[:-5] + 'qza')
        artifact = Artifact.import_data(
            'FeatureData[Taxonomy]', ref_taxa_fp,
            view_type='HeaderlessTSVTaxonomyFormat')
        artifact.save(ref_taxa_fp[:-3] + 'qza')


def test_cross_validated_sequences(data_dir):
    '''confirm that test (query) taxa IDs are not in training (ref) set, but
    that all taxonomy strings are.
    '''
    simulated_dir = join(data_dir, 'cross-validated')
    for db_iter_dir in glob(join(simulated_dir, '*-iter*')):
        query_taxa = import_to_list(
            join(db_iter_dir, 'query_taxa.tsv'), field=0)
        ref_taxa = import_to_list(
            join(db_iter_dir, 'ref_taxa.tsv'), field=0)
        common = set(query_taxa).intersection(set(ref_taxa))
        if common:
            print('sequences in training and test sets:')
            print(common)


def test_novel_taxa_datasets(data_dir):
    '''confirm that test (query) taxa IDs and taxonomies are not in training
    (ref) set, but sister branch taxa are.
    '''
    novel_dir = join(data_dir, 'novel-taxa-simulations')
    for db_iter_dir in glob(join(novel_dir, '*-L*-iter*')):
        query_taxa = import_taxonomy_to_dict(join(db_iter_dir,
                                                  'query_taxa.tsv'))
        ref_taxa = import_taxonomy_to_dict(join(db_iter_dir,
                                                'ref_taxa.tsv'))
        for key, value in query_taxa.items():
            if key in ref_taxa:
                print('key duplicate:', basename(db_iter_dir), key)
            if value in ref_taxa.values():
                print('value duplicate:', basename(db_iter_dir), value)


def recall_novel_taxa_dirs(data_dir, databases, iterations,
                           ref_seqs='ref_seqs.fasta', ref_taxa='ref_taxa.tsv',
                           max_level=6, min_level=0, multilevel=True):
    '''Given the number of iterations and database names, create list of
    directory names, and dict of reference seqs and reference taxa.
    data_dir = base directory containing results directories
    databases = names of ref databases used to generate novel taxa datasets
    iterations = number of iterations set during novel taxa dataset generation
    ref_seqs = filepath of reference sequences used for assignment
    ref_taxa = filepath of reference taxonomies used for assignment
    max_level = top level INDEX in RANGE to recall
    min_level = bottom level INDEX in RANGE to recall
                e.g., max_level=6, min_level=0 generates 1,2,3,4,5,6
    multilevel = whether taxa assignments should be iterated over multiple
                 taxonomic levels (as with novel taxa). Set as False if taxa
                 assignment should not be performed at multiple levels, e.g.,
                 for simulated community analysis. Levels must still be set to
                 iterate across a single level, e.g., max_level=6, min_level=5
    '''
    dataset_reference_combinations = list()
    reference_dbs = dict()
    for database in databases:
        for level in range(max_level, min_level, -1):
            for iteration in range(0, iterations):
                if multilevel is True:
                    dataset_name = '{0}-L{1}-iter{2}'.format(database,
                                                             level,
                                                             iteration)
                else:
                    dataset_name = '{0}-iter{1}'.format(database, iteration)
                dataset_reference_combinations.append((dataset_name,
                                                       dataset_name))
                reference_dbs[dataset_name] = (join(data_dir, dataset_name,
                                                    ref_seqs),
                                               join(data_dir, dataset_name,
                                                    ref_taxa))
    return dataset_reference_combinations, reference_dbs


# tag for removal — if match == recall, do we need to keep LCA?
def find_last_common_ancestor(taxon_a, taxon_b):
    '''Given two taxonomy strings, find the first level at which taxa mismatch.
    '''
    i = -1
    for i, (ta, tb) in enumerate(zip(taxon_a.split(';'), taxon_b.split(';'))):
        if ta != tb:
            return i
    return i + 1


def evaluate_classification(obs_taxon, exp_taxon):
    '''Given an observed and actual taxonomy string corresponding to a cross-
    validated simulated community, score as match, overclassification,
    underclassification, or misclassification'''

    if obs_taxon == exp_taxon:
        return 'match'
    if exp_taxon.startswith(obs_taxon) or \
            obs_taxon in ('Unclassified', 'Unassigned', 'No blast hit'):
        return 'underclassification'
    if obs_taxon.startswith(exp_taxon):
        return 'overclassification'
    return 'misclassification'


def load_taxa(obs_fp, level=slice(0, 7), field=1, sort=True):
    '''Mount observed/expected taxonomy observations.
    obs_fp: path
        Input file containing taxonomy strings and IDs.
    level: slice
        Taxonomic range of interest. 0 = kingdom, 1 = phylum, 2 = class,
        3 = order, 4 = family, 5 = genus, 6 = species.
        Must use slice notation. For species, use level=slice(6, 7)
    field: int
        ab-delimited field containing taxonomy strings.
    '''
    raw_obs = accept_list_or_file(obs_fp)
    if sort:
        raw_obs = sorted(raw_obs)
    obs = extract_taxa_names(raw_obs, field=field, level=level)
    obs = [';'.join([l.strip() for l in line.split(';')]) for line in obs
           if not line.startswith(('taxonomy', 'Taxon'))]
    return obs


def count_records(record_counter, record_name, line_count):
    '''Tally ratio of results in record counter.
    record_counter: counter
        Name of counter containing record_name and line_count counts.
    record_name: str
        Key name for record to tally.
    count: str
        Key name for record of total lines, used as denominator.
    '''
    count = record_counter[line_count]
    try:
        ratio = record_counter[record_name] / count
    except ZeroDivisionError:
        ratio = 0
    return ratio


def precision_recall_fscore(exp, obs, sample_weight=None, exclude=None):
    # precision, recall, fscore, calculated using microaveraging
    if exclude is None:
        exclude = []
    if sample_weight is None:
        sample_weight = [1]*len(exp)
    tp, fp, fn = 0, 0, 0
    for e, o, w in zip(exp, obs, sample_weight):
        if e in exclude:
            continue
        classification = evaluate_classification(o, e)
        if classification == 'match':
            # tp for the true class, the rest are tn
            tp += w
        elif classification == 'underclassification':
            # fn for the true class
            fn += w
            # no fp for the predicted class, because it was right to some level
            # the rest are tn
        else:
            # fp the the predicted class
            fp += w
            # fn for the true class, the rest are tn
            fn += w
    if tp == 0:
        return 0, 0, 0
    p = tp / (tp + fp)
    r = tp / (tp + fn)
    f = 2.*p*r / (p + r)
    return p, r, f


def compute_prf(exp, obs, test_type='cross-validated',
                l_range=range(1, 7), sample_weight=None, exclude=None):
    '''Compute precision, recall, and F-measure using sklearn.
    exp_taxa: list
        Expected observations for each sample (sequence).
    obs_taxa: list
        Predicted observations for each sample (sequence).
    avg: str
        Score averaging method using in sklearn. 'micro', 'weighted', or
        'macro'.
    test_type: str
        'novel-taxa' or 'cross-validated'.
    l_range: range
        Range of taxonomic levels to test is test_type = 'cross-validated'.
    sample_weight: array-like of shape = [n_samples], optional
        Sample weights.
    '''

    if test_type in ('mock', 'novel-taxa'):
        p, r, f = precision_recall_fscore(
            exp, obs, sample_weight=sample_weight, exclude=exclude)
    elif test_type == 'cross-validated':
        # initialize p/r/f as lists of 0s, representing each taxonomic level.
        p, r, f = [0] * 7, [0] * 7, [0] * 7
        # iterate over multiple taxonomic levels
        for lvl in l_range:
            _obs = extract_taxa_names(obs, level=slice(0, lvl+1))
            _exp = extract_taxa_names(exp, level=slice(0, lvl+1))
            p[lvl], r[lvl], f[lvl] = precision_recall_fscore(
                _exp, _obs, sample_weight=sample_weight, exclude=exclude)
    else:
        raise ValueError('test_type must be "novel-taxa" or "cross-validated" '
                         'or "mock".')

    return p, r, f


def load_prf(obs_fp, exp_fp, level=slice(0, 7), sort=True):
    exp_taxa = load_taxa(exp_fp, level=level, sort=sort)
    obs_taxa = load_taxa(obs_fp, level=level, sort=sort)

    # raise error if obs_taxa and exp_taxa are not same length
    if len(obs_taxa) != len(exp_taxa):
        raise RuntimeError(
            'Lengths of expected and observed taxa do not match. '
            'Check inputs: {0}, {1}'.format(obs_fp, exp_fp))

    return exp_taxa, obs_taxa


def novel_taxa_classification_evaluation(results_dirs, expected_results_dir,
                                         summary_fp, test_type='novel-taxa'):
    '''Input glob of novel taxa results, receive a summary of accuracy results.
    results_dirs = list or glob of novel taxa observed results in format:
                    precomputed_results_dir/dataset_id/method_id/params_id/
    expected_results_dir = directory containing expected novel-taxa results in
                    format:
                    expected_results_dir/dataset_id/method_id/params_id/
    summary_fp = filepath to contain summary of results
    test_type = 'novel-taxa' or 'cross-validated'

    Returns results as df, in addition to printing summary_fp
    '''
    results = []

    for results_dir in results_dirs:
        fields = results_dir.split(sep)
        dataset_id, method_id, params_id = fields[-3], fields[-2], fields[-1]

        if test_type == 'novel-taxa':
            index = dataset_id.split('-L')[0]
            level = int(dataset_id.split('-')[2].lstrip('L').strip())
            iteration = dataset_id.split('iter')[1]
        elif test_type == 'cross-validated':
            index, iteration = dataset_id.split('-iter')
            level = 6

        # import observed and expected taxonomies to list; order both by ID
        obs_fp = join(results_dir, 'query_tax_assignments.txt')
        exp_fp = join(expected_results_dir, dataset_id, 'query_taxa.tsv')
        exp_taxa, obs_taxa = load_prf(obs_fp, exp_fp)

        p, r, f = compute_prf(exp_taxa, obs_taxa, test_type=test_type)

        # Create empty list of levels at which first mismatch occurs
        mismatch_level_list = [0] * 8
        log = ['dataset\tlevel\titeration\tmethod\tparameters\
               \tobserved_taxonomy\texpected_taxonomy\tresult\tmismatch_level\
               \tPrecision\tRecall\tF-measure']

        # loop through observations, store results to counter
        record_counter = Counter()
        for obs, exp in zip(obs_taxa, exp_taxa):
            # Find shallowest level of mismatch
            mismatch_level = find_last_common_ancestor(obs, exp)
            mismatch_level_list[mismatch_level] += 1

            # evaluate novel taxa classification
            result = evaluate_classification(obs, exp)

            record_counter.update({'line_count': 1})
            record_counter.update({result: 1})
            log.append('\t'.join(map(str, [index, level, iteration,
                                           method_id, params_id,
                                           obs, exp, result,
                                           mismatch_level, p, r, f])))

        # Create log file
        log_fp = join(results_dir, 'classification_accuracy_log.tsv')
        export_list_to_file(log, log_fp)

        # tally score ratios
        match_ratio = count_records(record_counter, 'match', 'line_count')
        overclass = count_records(record_counter, 'overclassification',
                                  'line_count')
        underclass = count_records(record_counter, 'underclassification',
                                   'line_count')
        misclass = count_records(record_counter, 'misclassification',
                                 'line_count')

        # add everything to results
        results.append((index, level, iteration, method_id, params_id,
                        match_ratio, overclass, underclass, misclass,
                        mismatch_level_list, p, r, f))

    # send to dataframe, write to summary_fp
    result = pd.DataFrame(results, columns=["Dataset", "level", "iteration",
                                            "Method", "Parameters",
                                            "match_ratio",
                                            "overclassification_ratio",
                                            "underclassification_ratio",
                                            "misclassification_ratio",
                                            "mismatch_level_list", "Precision",
                                            "Recall", "F-measure"])
    result.to_csv(summary_fp)
    return result


# tag for modification: if we remove match/LCA, remove support from this fn
def extract_per_level_accuracy(df, columns=['Precision', 'Recall', 'F-measure',
                                            'mismatch_level_list']):
    '''Generate new pandas dataframe, containing match ratios for taxonomic
    assignments at each taxonomic level. Extracts mismatch_level_list from a
    dataframe and splits this into separate df entries for plotting.

    df: dataframe
        pandas dataframe
    column: list
        column names containing mismatch_level_list or other lists to be
        separated into multiple dataframe entries for plotting.

        mismatch_level_list reports mismatches at each level of taxonomic
        assignment (8 levels).

        Currently levels  are hardcoded, but could be adjusted
        below in lines:
            for level in range(1, 7):
    '''
    results = []

    for index, data in df.iterrows():
        for level in range(1, 7):
            level_results = []
            col_names = []
            for column in columns:
                # If using precomputed results, mismatch_level_list is imported
                # as string, hence must be converted back to list of integers.
                if isinstance(data[column], str):
                    col = list(map(float, data[column].strip('[]').split(',')))
                else:
                    col = data[column]
                # 'mismatch_level_list' contains level of first mismatch for
                # each observation; hence, matches at level L = total
                # observations - cumulative mismatches / total observations
                if column == 'mismatch_level_list':
                    linecount = sum(col)
                    col_names.append("match_ratio")
                    cumulative_mismatches = sum(col[0:level+1])
                    if cumulative_mismatches < linecount:
                        score = (linecount - cumulative_mismatches) / linecount
                    else:
                        score = col[0]
                # Otherwise just extract score at level index.
                else:
                    score = col[level]
                    col_names.append(column)

                # add column score for level to level_results
                level_results.append(score)

            results.append((data['Dataset'], level, data['iteration'],
                            data['Method'], data['Parameters'],
                            *[s for s in level_results]))

    result = pd.DataFrame(results, columns=["Dataset", "level", "iteration",
                                            "Method", "Parameters",
                                            *[s for s in col_names]])
    return result


def runtime_make_test_data(seqs_in, results_dir, sampling_depths):
    '''Repeatedly subsample a fasta sequence file at multiple sequence depths
    to generate query/test data for testing method runtimes.

    seqs_in: path
        fasta format reference sequences.
    results_dir: path
        Output directory.
    sampling_depths: list of integers
        Number of sequences to subsample from seqs.
    '''
    if not exists(results_dir):
        makedirs(results_dir)

    seqs = [seq for seq in io.read(seqs_in, format='fasta')]
    for depth in sampling_depths:
        subset = sample(seqs, depth)
        tmpfile = join(results_dir, str(depth)) + '.fna'
        with open(tmpfile, "w") as output_fasta:
            for s in subset:
                s.write(output_fasta, format='fasta')


def runtime_make_commands(input_dir, results_dir, methods,
                          ref_taxa, sampling_depths, num_iters=1,
                          subsample_ref=True):
    '''Generate list of commands to benchmark method runtimes.

    input_dir: path
        Input directory, containing query/ref sequences.
    results_dir: path
        Output directory.
    methods: dict
        Dictionary of method:parameters pairs in format:
            {'method' : (command-template, method-specific-parameters)}
    ref_taxa: path
        Taxonomy map for ref sequences in tab-separated format:
            seqID   ACGTGTAGTCGATGCTAGCTACG
    sampling_depths: list of integers
        Number of sequences to subsample from seqs.
    num_iters: int
        Number of iterations to perform.
    subsample_ref: bool
        If True (default), ref seqs are subsampled at depths defined in
        sampling_depths, and query seqs default to smallest depth. If false,
        query seqs are subsampled at these depths, and ref defaults to largest
        sampling depth.
    '''

    commands = []
    for iteration in range(num_iters):
        for method, template in methods.items():
            for depth in sampling_depths:
                # default: subsample ref seqs, query = smallest sample
                if subsample_ref is True:
                    q_depth = str(sampling_depths[0])
                    r_depth = str(depth)
                # or subsample query seqs, ref = largest sample
                else:
                    q_depth = str(depth)
                    r_depth = str(sampling_depths[-1])
                query = join(input_dir, q_depth) + '.fna'
                ref = join(input_dir, r_depth) + '.fna'
                command = (template[0].format(results_dir, query, ref,
                                              ref_taxa, method, template[1]),
                           method, q_depth, r_depth, iteration)
                commands.append(command)
    return commands


def clock_runtime(command, results_fp, force=False):
    '''Execute a command and record the runtime.

    command: str
        Command to be executed.
    results_fp: path
        Output file
    force: bool
        Overwrite results? If false, will append to any existing results
    '''
    if force is True:
        remove(results_fp)

    _command, method, q_frac, r_frac, iteration = command
    start = time()
    system(_command)
    end = time()
    results = [method, q_frac, r_frac, iteration, end - start]
    with open(results_fp, 'a') as timeout:
        timeout.write('\t'.join(map(str, results)) + '\n')
