#!/usr/bin/env python


# ----------------------------------------------------------------------------
# Copyright (c) 2016--, tax-credit development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


from tax_credit.taxa_manipulator import (compile_reference_taxa,
                                         filter_sequences,
                                         convert_tsv_to_biom,
                                         extract_fasta_ids,
                                         string_search)
from re import search
from random import choice, random
from shutil import copyfile
from os import makedirs
from os.path import join, exists
import pandas as pd


def generate_simulated_communities(sim_dir, datasets, reference_dbs,
                                   strain_max=5):
    '''Generate simulated communities based on natural community composition
    and reference sequences matching the natural taxonomic assignments.
    __________
    Parameters
    __________
    sim_dir: path
        directory containing simulated community dataset directories.
    datasets: list
        Tuples containg directory name and name of reference db that is a key
        in reference_dbs: (directory name, reference name)
        The combination of "<directory name>/<reference name>" indicated in the
        tuple values must include 'expected-composition.txt', containing the
        expected composition of each sample in the format:
        #SampleID   taxa1   taxa2   taxa3
        sample1     0.1     0.3     0.6
        sample2     0.4     0.3     0.3

    reference_dbs: dict
        Contains paths to reference seqs (fasta) and taxonomy (tab-separated
        format with sequence ID <tab> semicolon-separated taxonomy annotation
        on each line).
        {name: (sequences, taxonomy)}
    strain_max: int
        Each taxonomy in expected-composition will be matched against up to
        strain_max sequence IDs in the reference taxonomy, and the
        corresponding sequences will be used as representative sequences for
        that taxonomy.
    __________
    Returns
    __________
    Writes to file:
    1. Simulated compositions in tsv format as for expected-composition.
    2. biom tables for simulated-compositions and expected-composition.
    3. A fasta file containing "representative sequences" for each "OTU",
        i.e., reference sequences matching the expected taxa, up to strain_max
        sequences per taxonomy.

    '''
    for dataset, ref_data in datasets:
        # locate source data
        community_dir = join(sim_dir, dataset)

        # extract taxa names from expected-composition
        exp_comp_fp = join(community_dir, 'expected-composition.txt')
        comp = pd.DataFrame.from_csv(exp_comp_fp, sep='\t')
        # extract taxa and per-sample frequencies
        taxa = {t: {s: comp.loc[s, t]
                    for s in comp.index} for t in comp.columns}

        # pull matching reference seqs
        ref_seqs_fp, ref_taxa_fp = reference_dbs[ref_data]
        # First, make sure taxa represent only clean read seqs
        valid_t = '^' + '$|^'.join(list(extract_fasta_ids(ref_seqs_fp))) + '$'
        read_taxa = string_search(ref_taxa_fp, valid_t, discard=False,
                                  field=slice(0, 1), delim='\t')
        ref_taxa = compile_reference_taxa(read_taxa)

        rep_seq_ids = []
        sim_comp = {}
        c_match = 0
        c_close = 0
        for t, vals in taxa.items():

            # compile matching ref seqs
            taxa_matches = []
            # attempt exact match
            if t in ref_taxa:
                taxa_matches.append((ref_taxa[t]))
                c_match += 1
            # otherwise pull all partial matches
            else:
                taxa_matches = [ref_taxa[rt]
                                for rt in ref_taxa.keys()
                                if search(t, rt)]
                c_close += 1

            # randomly sample taxa_matches for seq ids twice:
            # first for taxa, then sample taxa choice for id.
            try:
                taxa_ids = list({choice(choice(taxa_matches))
                                 for n in range(strain_max)})
            except IndexError:
                raise IndexError("one or more taxa do not match reference.")
            rep_seq_ids = rep_seq_ids + taxa_ids

            # generate simulated composition
            for t_id in taxa_ids:
                # randomly split val frequency into each component strain
                sim_comp[t_id] = {}
                for s_id, val in vals.items():
                    rands = [random() for n in range(strain_max)]
                    rands = [n / sum(rands) for n in rands]
                    sim_comp[t_id][s_id] = val * rands.pop()

        # create dataframe of simulated composition
        simulation = pd.DataFrame.from_dict(sim_comp)
        sim_comp_fp = join(community_dir, 'simulated-composition.txt')
        sim_biom_fp = join(community_dir, 'simulated-composition.biom')
        exp_biom_fp = join(community_dir, 'expected-composition.biom')
        simulation.to_csv(sim_comp_fp, sep='\t', index_label='#SampleID')
        convert_tsv_to_biom(sim_comp_fp, sim_biom_fp)
        convert_tsv_to_biom(exp_comp_fp, exp_biom_fp, obs_ids_to_metadata=True)

        # generate rep seqs
        filter_sequences(
            ref_seqs_fp, join(community_dir, 'simulated-seqs.fna'),
            rep_seq_ids, keep=True)

        print("{0}: {1} matches and {2} near matches.".format(
            dataset, c_match, c_close))


def copy_expected_composition(source_dir, dataset_reference_combinations,
                              dest_dir, fn='expected-composition.biom'):
    '''Copy expected composition files to precomputed results.'''
    for dataset, ref in dataset_reference_combinations:
        source_path = join(source_dir, dataset, fn)
        dest_path = join(dest_dir, dataset, ref, 'expected')
        if not exists(dest_path):
            makedirs(dest_path)
        copyfile(source_path, join(dest_path, fn))
