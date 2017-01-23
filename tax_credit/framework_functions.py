#!/usr/bin/env python


# ----------------------------------------------------------------------------
# Copyright (c) 2016--, tax-credit development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from sys import exit
from os.path import join, exists, split, sep, expandvars, basename, splitext
from os import makedirs, remove, system
from glob import glob
from itertools import product
from shutil import rmtree, move
from random import sample
from biom import load_table
from biom.cli.util import write_biom_table
from biom.parse import MetadataMap
from skbio.alignment import local_pairwise_align_ssw
from skbio import io, DNA
from time import time
import pandas as pd
from collections import Counter
from sklearn.metrics import precision_recall_fscore_support
from tax_credit.taxa_manipulator import (accept_list_or_file,
                                         import_taxonomy_to_dict,
                                         export_list_to_file,
                                         extract_rownames,
                                         filter_sequences,
                                         extract_fasta_ids,
                                         string_search,
                                         trim_taxonomy_strings,
                                         unique_lines,
                                         branching_taxa,
                                         stratify_taxonomy_subsets,
                                         extract_taxa_names)


def parameter_sweep(data_dir, results_dir, reference_dbs,
                    dataset_reference_combinations,
                    method_parameters_combinations,
                    command_template=None, infile='rep_set.fna',
                    output_name='rep_set_tax_assignments.txt', force=False):
    '''Create list of commands from input dictionaries of method_parameter and
    dataset_reference_combinations
    '''

    if command_template is None:
        # default to qiime1 command template
        command_template = "mkdir -p {0} ; assign_taxonomy.py -v -i {1} -o {0}\
                            -r {2} -t {3} -m {4} {5} --rdp_max_memory 16000"

    commands = []
    for dataset, reference in dataset_reference_combinations:
        query_seqs = join(data_dir, dataset, infile)
        reference_seqs, reference_tax = reference_dbs[reference]
        for method, parameters in method_parameters_combinations.items():
            parameter_ids = sorted(parameters.keys())
            for parameter_combination in product(*[parameters[id_]\
                                                 for id_ in parameter_ids]):
                parameter_comb_id = ':'.join(map(str, parameter_combination))
                parameter_output_dir = join(results_dir, dataset, reference,
                                            method, parameter_comb_id)
                if not exists(join(parameter_output_dir, output_name))\
                  or force is True:
                    parameter_str = ' '.join(['--%s %s' % e for e in zip(\
                                        parameter_ids, parameter_combination)])
                    command = command_template.format(parameter_output_dir,
                                                      query_seqs,
                                                      reference_seqs,
                                                      reference_tax, method,
                                                      parameter_str)
                    commands.append(command)
    return commands


def add_metadata_to_biom_table(biom_input_fp, taxonomy_map_fp, biom_output_fp):
    '''Load biom, add metadata, write to new table'''
    newbiom = load_table(biom_input_fp)
    metadata = MetadataMap.from_file(taxonomy_map_fp,
                                     header=['Sample ID', 'taxonomy', 'c'])
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
                   junk='__;|__$|_sp_|unknown|unidentified'):

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

    clean_taxa = join(db_dir,
                      '{0}_clean.tsv'.format(basename(splitext(taxa_in)[0])))
    clean_fasta = join(db_dir,
                      '{0}_clean.fasta'.format(basename(splitext(seqs_in)[0])))

    # Remove empty taxa from ref taxonomy
    taxa = string_search(taxa_in, junk, discard=True)
    # Remove brackets (and other special characters causing problems)
    clean_list = [line.translate(str.maketrans('', '', '[]()'))\
                  for line in taxa]
    export_list_to_file(clean_list, clean_taxa)
    # Remove corresponding seqs from ref fasta
    filter_sequences(seqs_in, clean_fasta, clean_taxa, keep=True)

    return clean_taxa, clean_fasta


def find_primer_site(query_seq, target_seq):
    '''Align primer to target DNA using striped Smith-Waterman, return site'''
    # Increase gap extend penalty to avoid large gaps, short alignments
    # try:
    aln, score, start_end_positions = local_pairwise_align_ssw(query_seq,
                                        target_seq, gap_extend_penalty=10)
    # except IndexError:
    #    continue
    return start_end_positions


def extract_amplicons(seqs_in, amplicons_out, reads_out, fwd_primer,
                      rev_primer, read_length, min_read_length=80):
    '''Generate artificial amplicons and reads from a set of fasta sequences'''

    with open(amplicons_out, 'w') as amplicons:
        with open(reads_out, 'w') as reads:
            for seq in io.read(seqs_in, format='fasta'):
                # Align forward and reverse primers onto input sequence
                seq_DNA = DNA(seq, lowercase=True)
                fwd_primer_start = find_primer_site(fwd_primer, seq_DNA)[1][0]
                rev_primer_end = find_primer_site(rev_primer, seq_DNA)[1][1]
                # Trim seq to amplicon
                amplicon = seq[fwd_primer_start:rev_primer_end]
                if len(amplicon) != 0:
                    amplicon.write(amplicons)
                # Trim amplicon to sequence read length
                read = amplicon[0:read_length]
                if len(read) >= min_read_length:
                    read.write(reads)


def seq_count(infile):
    '''Count sequences in fasta file'''
    with open(infile, "r") as f:
        count = sum(1 for line in f) / 2
    return count


def generate_simulated_datasets(dataframe, data_dir, read_length, iterations):
    '''From a dataframe of sequence reference databases, build training/test
    sets of "novel taxa" queries, taxonomies, and reference seqs/taxonomies.

    dataframe = pandas dataframe
    '''
    for index, data in dataframe.iterrows():
        db_dir = join(data_dir, 'ref_dbs', data['Reference id'])
        if not exists(db_dir):
            makedirs(db_dir)

        # Clean taxonomy/sequences, to remove empty/ambiguous taxonomies
        clean_fasta = join(db_dir, '{0}_clean.fasta'.format(basename(splitext(\
                                             data['Reference file path'])[0])))
        if not exists(clean_fasta):
            clean_taxa, clean_fasta = clean_database(data['Reference tax path'],
                                                    data['Reference file path'],
                                                    db_dir)

        # Trim reference sequences to amplicon target
        primer_pair = '{0}-{1}'.format(data['Fwd primer id'],
                                       data['Rev primer id'])
        base, ext = splitext(clean_fasta)
        amplicons_fp = join(db_dir, '{0}_{1}{2}'.format(base, primer_pair, ext))
        simulated_reads_fp = join(db_dir, '{0}_{1}_trim{2}{3}'.format(base,
                                                               primer_pair,
                                                               read_length,
                                                               ext))

        #amplicons_fp = join(db_dir, "simulated_amplicons.fna")
        #simulated_reads_fp = join(db_dir, "simulated_reads.fna")
        if not exists(simulated_reads_fp):
            extract_amplicons(clean_fasta, amplicons_fp, simulated_reads_fp,
                          DNA(data['Fwd primer']), DNA(data['Rev primer']),
                          read_length, min_read_length=80)
        else:
            print('simulated reads and amplicons exist: skipping extraction')

        # Filter taxonomy strings to match read sequences
        valid_taxa = '^' + '$|^'.join(list(extract_fasta_ids(\
                                            simulated_reads_fp))) + '$'
        read_taxa = string_search(clean_taxa, valid_taxa, discard = False,
                                  field=slice(0, 1), delim='\t')

        # Print sequence counts to see how many seqs were filtered
        print(index, 'Sequence Counts')
        print('Raw Fasta:           ', seq_count(data['Reference file path']))
        print('Clean Fasta:         ', seq_count(clean_fasta))
        print('Simulated Amplicons: ', seq_count(amplicons_fp))
        print('Simulated Reads:     ', seq_count(simulated_reads_fp))

        # Generate novel query and reference seqs/taxa pairs
        novel_dir = join(data_dir, 'novel-taxa-simulations')
        generate_novel_sequence_sets(read_taxa, simulated_reads_fp, index,
                                     iterations, novel_dir)

        # Generate simulated community query and reference seqs/taxa pairs
        simulated_dir = join(data_dir, 'simulated-community')
        generate_simulated_communities(read_taxa, simulated_reads_fp, index,
                                     iterations, simulated_dir)


def generate_novel_sequence_sets(read_taxa, simulated_reads_fp, index,
                                 iterations, data_dir):
    '''Generate paired query/reference fastas and taxonomies for novel taxa
    analysis, given an input of simulated amplicon taxonomies (read_taxa)
    and fastas (simulated_reads_fp), the index (database) name, # of
    iterations to perform (cross-validated data subsets), the output dir
    (data_dir).
    read_taxa = list or file of taxonomies corresponding to simulated_reads_fp
    simulated_reads_fp = simulated amplicon reads (fasta format file)
    index = reference database name
    iterations = number of subsets to create
    data_dir = base output directory to contain simulated datasets
    '''
    # Remove non-branching taxa
    for level in range(6, 0, -1):
        # Trim taxonomy strings to level X
        expected_taxa = trim_taxonomy_strings(read_taxa, level)
        # sample unique representative taxa
        unique_taxa = unique_lines(expected_taxa, mode='n', field=1)
        # Remove non-branched taxa
        branched_taxa = branching_taxa(unique_taxa, field=level)

        # Count unique and braching taxa, just for fun
        print('{0} level {1} contains {2} unique and {3} branched taxa\
              '.format(index, level, len(unique_taxa), len(branched_taxa)))

    # Create QUERY / REF pairs
        # Create QUERY TAXONOMY as subsets of branched taxa lists
        basename = '{0}-L{1}'.format(index, level)
        # stratify at level-1, since we want to stratify the branches, not tips
        stratify_taxonomy_subsets(branched_taxa, iterations, data_dir,
                                  basename, level=level-1)

        # Create REF and TAXONOMY files paired to QUERY files:
        for iteration in range(0, iterations):
            db_iter_dir = join(data_dir, '{0}-iter{1}'.format(basename,
                                                              iteration))
            novel_query_taxonomy_fp = join(db_iter_dir,
                                                'query_taxa.tsv')
            novel_query_fp = join(db_iter_dir, 'query.fasta')
            novel_ref_fp = join(db_iter_dir, 'ref_seqs.fasta')
            novel_ref_taxonomy_fp = join(db_iter_dir, 'ref_taxa.tsv')

            # 1) Create REF TAXONOMY from list of taxonomies that
            #    DO NOT match QUERY taxonomies
            name_list = '^' + '$|^'.join(list(set(extract_taxa_names(\
                         novel_query_taxonomy_fp, slice(1, level+1))))) + '$'
            novel_ref_taxonomy = string_search(read_taxa, name_list,
                                             discard = True,
                                             field = slice(1, level+1))
            export_list_to_file(novel_ref_taxonomy, novel_ref_taxonomy_fp)

            # 2) Create REF: Filter ref database to contain only seqs that
            #    match non-matching taxonomy strings
            filter_sequences(simulated_reads_fp, novel_ref_fp,
                             novel_ref_taxonomy_fp, keep=True)

            # 3) Create QUERY: Filter ref database to contain only seqs
            #    that match QUERY TAXONOMY
            filter_sequences(simulated_reads_fp, novel_query_fp,
                             novel_query_taxonomy_fp, keep=True)


def generate_simulated_communities(read_taxa, simulated_reads_fp, index,
                                   iterations, data_dir):
    '''Generates simulated community files (fasta and taxonomy) as subsets of
    simulated amplicons/taxa for cross-validated taxonomy assignment. Selects
    duplicated taxa names, evenly allocates these among subsets as query taxa
    (test set), generates ref taxa (training set) that do not match query fasta
    IDs, and creates fasta files to match each of these sets.
    read_taxa = list or file of taxonomies corresponding to simulated_reads_fp
    simulated_reads_fp = simulated amplicon reads (fasta format file)
    index = reference database name
    iterations = number of subsets to create
    data_dir = base output directory to contain simulated datasets
    '''
    # Subset amplicons so that taxa are evenly distributed between train/test
    duplicated_taxa = unique_lines(read_taxa, mode='d', field=1)
    stratify_taxonomy_subsets(duplicated_taxa, iterations, data_dir, index,
                              level=1, delim='\t')
    # generate pairs of train/test fastas/taxonomies for each CV subset
    for iteration in range(0, iterations):
        db_iter_dir = join(data_dir, '{0}-iter{1}'.format(index, iteration))
        query_taxa_fp = join(db_iter_dir, 'query_taxa.tsv')
        query_fp = join(db_iter_dir, 'query.fasta')
        ref_fp = join(db_iter_dir, 'ref_seqs.fasta')
        ref_taxa_fp = join(db_iter_dir, 'ref_taxa.tsv')

        # 1) Create REF taxa that do not match query IDs
        ids = '^' + '$|^'.join(list(extract_rownames(query_taxa_fp))) + '$'
        ref_taxa = string_search(read_taxa, ids, discard=True,
                                 field=slice(0,1), delim='\t')
        export_list_to_file(ref_taxa, ref_taxa_fp)
        # 2) Create REF: seqs that match ref taxonomy
        filter_sequences(simulated_reads_fp, ref_fp, ref_taxa_fp, keep=True)
        # 3) Create QUERY: seqs that match QUERY TAXONOMY
        filter_sequences(simulated_reads_fp, query_fp, query_taxa_fp, keep=True)


def test_simulated_communities(dataframe, data_dir, iterations):
    '''confirm that test (query) taxa IDs are not in training (ref) set, but
    that all taxonomy strings are.
    '''
    simulated_dir = join(data_dir, 'simulated-community')
    for index, data in dataframe.iterrows():
        for iteration in range(0, iterations):
            db_iter_dir = join(simulated_dir, '{0}-iter{1}'.format(index,
                                                                   iteration))
            query_taxa = import_taxonomy_to_dict(join(db_iter_dir,
                                                           'query_taxa.tsv'))
            ref_taxa = import_taxonomy_to_dict(join(db_iter_dir,
                                                         'ref_taxa.tsv'))
            for key, value in query_taxa.items():
                if key in ref_taxa:
                    print('key duplicate: ', key)
                if value not in ref_taxa.values():
                    print('missing value: ', value)


def test_novel_taxa_datasets(dataframe, data_dir, iterations):
    '''confirm that test (query) taxa IDs and taxonomies are not in training
    (ref) set, but sister branch taxa are.
    '''
    novel_dir = join(data_dir, 'novel-taxa-simulations')
    for index, data in dataframe.iterrows():
        for level in range(6, 0, -1):
            for iteration in range(0, iterations):
                db_iter_dir = join(novel_dir, '{0}-L{1}-iter{2}'.format(index,
                                                                    level,
                                                                    iteration))
                query_taxa = import_taxonomy_to_dict(join(db_iter_dir,
                                                             'query_taxa.tsv'))
                ref_taxa = import_taxonomy_to_dict(join(db_iter_dir,
                                                             'ref_taxa.tsv'))
                taxa = [t.split(';')[level-1] for t in ref_taxa.values()]
                for key, value in query_taxa.items():
                    if key in ref_taxa:
                        print('key duplicate:', index, level, iteration, key)
                    if value in ref_taxa.values():
                        print('value duplicate:', index, level, iteration,
                              value)
                    if value.split(';')[level-1] not in taxa:
                        print('missing branch:', index, level, iteration,
                              value.split(';')[level-1])


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
                                                   ref_taxa)
                                              )
    return dataset_reference_combinations, reference_dbs


def find_last_common_ancestor(taxon_a, taxon_b):
    '''Given two taxonomy strings (input as lists of taxa names separated by
    level), find the first level at which taxa mismatch.
    '''
    mismatch_level = 0
    for taxa_comparison in zip(taxon_a, taxon_b):
        if taxa_comparison[0].strip() == taxa_comparison[1].strip():
            mismatch_level += 1
        else:
            break
    return mismatch_level


def evaluate_novel_taxa_classification(obs_taxa, exp_taxa, level):
    '''Given an observed and actual taxonomy string corresponding to a "novel"
    taxon, score as match, overclassification, underclassification, or
    misclassification'''
    # if observed has same assignment depth as expected and top level matches,
    # ==match. len(exp_taxa) - 1 because exp_taxa is actual taxonomy string,
    # L-1 is the actual expected taxonomy string
    if len(obs_taxa) == len(exp_taxa) - 1 \
    and obs_taxa[level - 1].strip() == exp_taxa[level - 1].strip():
        result = 'match'
    # if deeper and assignemnt at L-1 is correct, count as overclassification
    elif len(obs_taxa) >= len(exp_taxa) \
    and obs_taxa[level - 1].strip() == exp_taxa[level - 1].strip():
        result = 'overclassification'
    # if shallower and top-level assign correct, count as underclassification
    elif len(obs_taxa) < len(exp_taxa) - 1 \
    and obs_taxa[len(obs_taxa)-1].strip() == exp_taxa[len(obs_taxa)-1].strip()\
    or obs_taxa[0] == 'Unclassified' or obs_taxa[0] == 'Unassigned':
        result = 'underclassification'
    # Otherwise, count as misclassification
    else:
        result = 'misclassification'
    return result


def evaluate_cross_validated_classification(obs_taxa, exp_taxa):
    '''Given an observed and actual taxonomy string corresponding to a cross-
    validated simulated community, score as match, overclassification,
    underclassification, or misclassification'''
    # if  observed = expected, match
    if len(obs_taxa) == len(exp_taxa)\
        and obs_taxa[len(obs_taxa)-1].strip() == exp_taxa[len(obs_taxa)-1].strip():
            result = 'match'
    # if shallower and top-level assign correct, count as underclassification
    elif len(obs_taxa) < len(exp_taxa) \
    and obs_taxa[len(obs_taxa)-1].strip() == \
    exp_taxa[len(obs_taxa)-1].strip()\
    or obs_taxa[0] == 'Unclassified' or obs_taxa[0] == 'Unassigned':
        result = 'underclassification'
    # Otherwise, count as misclassification
    else:
        result = 'misclassification'
    return result


def load_taxa(obs_fp, level=slice(0, 7), field=1):
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
    obs = extract_taxa_names(sorted(accept_list_or_file(obs_fp)), field=field,
                             level=level)
    obs = [';'.join([l.strip() for l in line.split(';')]) for line in obs
           if not line.startswith('taxonomy')]
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


def compute_prf(exp, obs, avg='micro', test_type='cross-validated',
                l_range=range(1,7), level=6):
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
    level:
        Level of taxonomic assignment used if test_type = 'novel-taxa'
    '''

    prf = precision_recall_fscore_support

    # p = tp / (tp + fp)
    # r = tp / (tp + fn)
    # f = (2 * p * r) / (p + r)
    # with labels, null classifications can be FN but not TP or FP.
    # Hence, nulls affect recall but not precision.

    def lab(exp, obs, level):
        # use set or else multiple counts weight observations.
        # Remove all labels < level. Hence, underclassifcation becomes null,
        # and can only be FN, but overclassification is still counted as FP.
        return [t for t in set(exp + obs) if not len(t.split(';')) < level]

    if test_type == 'novel-taxa':
        exp = extract_taxa_names(exp, level=slice(0, level))
            # slice at level, so that exp=level-1 (otherwise exp = true label,
            # not novel taxa label)
        p, r, f, s = prf(exp, obs, average=avg, labels=lab(exp, obs, level))
    elif test_type == 'cross-validated':
        # initialize p/r/f as lists of 0s, representing each taxonomic level.
        p, r, f = [0] * 7, [0] * 7, [0] * 7
        # iterate over multiple taxonomic levels
        for level in l_range:
            _obs = extract_taxa_names(obs, level=slice(0, level+1))
            _exp = extract_taxa_names(exp, level=slice(0, level+1))
            # Here use level+1 to slice actual level
            p[level], r[level], f[level], s = prf(_exp, _obs, average='micro',
                                                  labels=lab(_exp, _obs,
                                                  level+1))
    else:
        print('FAIL: test_type must == "novel-taxa" or "cross-validated"')

    return p, r, f


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
        exp_taxa = load_taxa(exp_fp, level=slice(0, 7))
        obs_taxa = load_taxa(obs_fp, level=slice(0, 7))

        # Exit if obs_taxa and exp_taxa are not same length
        if len(obs_taxa) != len(exp_taxa):
            exit('''FAIL: Lengths of expected and observed taxa do not match.
                 Check inputs: {0}, {1}'''.format(obs_fp, exp_fp))

        p, r, f = compute_prf(exp_taxa, obs_taxa, test_type=test_type,
                              level=level)

        # Create empty list of levels at which first mismatch occurs
        mismatch_level_list = [0] * 8
        log = ['dataset\tlevel\titeration\tmethod\tparameters\
               \tobserved_taxonomy\texpected_taxonomy\tresult\tmismatch_level\
               \tPrecision\tRecall\tF-measure']

        # loop through observations, store results to counter
        record_counter = Counter()
        for obs, exp in zip(obs_taxa, exp_taxa):
            # access "expected taxonomy"
            o = obs.split(';')
            e = exp.split(';')

            # Find shallowest level of mismatch
            mismatch_level = find_last_common_ancestor(o, e)
            mismatch_level_list[mismatch_level] += 1

            # evaluate novel taxa classification
            if test_type == 'novel-taxa':
                result = evaluate_novel_taxa_classification(o, e, level)

            elif test_type == 'cross-validated':
                result = evaluate_cross_validated_classification(o, e)

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


def novel_taxa_classification_evaluation_old(results_dirs, expected_results_dir,
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

        obs_taxa_fp = join(results_dir, 'query_tax_assignments.txt')
        exp_taxa_fp = join(expected_results_dir, dataset_id, 'query_taxa.tsv')

        # Create empty list of levels at which first mismatch occurs
        mismatch_level_list = [0, 0, 0, 0, 0, 0, 0, 0]

        # import expected taxonomies to list
        expectations = import_taxonomy_to_dict(exp_taxa_fp)

        log = ['dataset\tlevel\titeration\tmethod\tparameters\tseq_id\
               \tobserved_taxonomy\texpected_taxonomy\tresult\tmismatch_level']

        # loop through observations, store results to counter
        record_counter = Counter()
        with open(obs_taxa_fp, 'r') as observations:
            for line in observations:
                if line.startswith("#"):
                    continue

                obs_id = line.split('\t')[0]
                obs_taxonomy = line.split('\t')[1]

                # access "expected taxonomy"
                obs_taxa = obs_taxonomy.split(';')
                exp_taxonomy = expectations[obs_id]
                exp_taxa = exp_taxonomy.split(';')

                # Find shallowest level of mismatch
                mismatch_level = find_last_common_ancestor(obs_taxa, exp_taxa)
                mismatch_level_list[mismatch_level] += 1

                # evaluate novel taxa classification
                if test_type == 'novel-taxa':
                    result = evaluate_novel_taxa_classification(obs_taxa,
                                                                exp_taxa,
                                                                level)

                elif test_type == 'cross-validated':
                    result = evaluate_cross_validated_classification(obs_taxa,
                                                                     exp_taxa)

                record_counter.update({'line_count': 1})
                record_counter.update({result: 1})
                log.append('\t'.join(map(str, [index, level, iteration,
                                               method_id, params_id, obs_id,
                                               obs_taxonomy,
                                               exp_taxonomy, result,
                                               mismatch_level])))

        # Create log file
        log_fp = join(results_dir, 'classification_accuracy_log.tsv')
        export_list_to_file(log, log_fp)

        # tally score ratios
        count = record_counter['line_count']
        try:
            match_ratio = record_counter['match'] / count
        except ZeroDivisionError:
            match_ratio = 0
        try:
            overclassification_ratio = record_counter['overclassification'] / count
        except ZeroDivisionError:
            overclassification_ratio = 0
        try:
            underclassification_ratio = \
                record_counter['underclassification'] / count
        except ZeroDivisionError:
            underclassification_ratio = 0
        try:
            misclassification_ratio = record_counter['misclassification'] / count
        except ZeroDivisionError:
            misclassification_ratio = 0

        results.append((index, level, iteration, method_id,
                        params_id, match_ratio,
                        overclassification_ratio, underclassification_ratio,
                        misclassification_ratio, mismatch_level_list))

    # send to dataframe, write to summary_fp
    result = pd.DataFrame(results, columns=["Dataset", "level", "iteration",
                                            "Method", "Parameters",
                                            "match_ratio",
                                            "overclassification_ratio",
                                            "underclassification_ratio",
                                            "misclassification_ratio",
                                            "mismatch_level_list"])
    result.to_csv(summary_fp)
    return result


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
                # If using precomputed results, mismatch_level_list is imported as
                # string, hence must be converted back to list of integers.
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


def clock_runtime(command, results_fp, force=True):
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
        timeout.write('\t'.join(map(str,results)) + '\n')
