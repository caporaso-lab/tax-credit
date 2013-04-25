#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout", "Zack Ellett", "Kyle Patnode"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"
__status__ = "Development"

"""Contains functions used in the multiple_assign_taxonomy.py script."""

import sys
from time import time
from os.path import basename, isdir, join, normpath, split, splitext
from shutil import rmtree
from qiime.util import add_filename_suffix, create_dir
from qiime.workflow.util import (call_commands_serially, generate_log_fp,
                                 no_status_updates, print_commands,
                                 print_to_stdout, WorkflowError,
                                 WorkflowLogger)

def assign_taxonomy_multiple_times(input_dirs, output_dir, assignment_methods,
        reference_seqs_fp, id_to_taxonomy_fp,
        confidences=None, e_values=None, rtax_modes=None,
        input_fasta_filename='rep_set.fna',
        clean_otu_table_filename='otu_table_mc2.biom',
        read_1_seqs_filename='seqs1.fna', read_2_seqs_filename='seqs2.fna',
        rdp_max_memory=4000, command_handler=call_commands_serially,
        status_update_callback=no_status_updates, force=False):
    """ Performs sanity checks on passed arguments and directories. Builds 
        commands for each method and sends them off to be executed. """
    ## Check if output directory exists
    try:
        create_dir(output_dir, fail_on_exist=not force)
    except OSError:
        raise WorkflowError("Output directory '%s' already exists. Please "
                "choose a different directory, or force overwrite with -f."
                % output_dir)

    logger = WorkflowLogger(generate_log_fp(output_dir))

    for input_dir in input_dirs:
        ## Make sure the input dataset directory exists.
        if not isdir(input_dir):
            raise WorkflowError("The input dataset directory '%s' does not "
                                "exist." % input_dir)

        input_dir_name = split(normpath(input_dir))[1]
        output_dataset_dir = join(output_dir, input_dir_name)
        input_fasta_fp = join(input_dir, input_fasta_filename)
        clean_otu_table_fp = join(input_dir, clean_otu_table_filename)
        read_1_seqs_fp = join(input_dir, read_1_seqs_filename)
        read_2_seqs_fp = join(input_dir, read_2_seqs_filename)

        logger.write("\nCreating output subdirectory '%s' if it doesn't "
                     "already exist.\n" % output_dataset_dir)
        create_dir(output_dataset_dir)

        for method in assignment_methods:
            ## Method is RDP
            if method == 'rdp':
                ## Check for execution parameters required by RDP method
                if confidences is None:
                    raise WorkflowError("You must specify at least one "
                                        "confidence level.")
                ## Generate command for RDP
                commands = _generate_rdp_commands(output_dataset_dir,
                                                  input_fasta_fp,
                                                  reference_seqs_fp,
                                                  id_to_taxonomy_fp,
                                                  clean_otu_table_fp,
                                                  confidences,
                                                  rdp_max_memory=rdp_max_memory)
                        
            ## Method is BLAST
            elif method == 'blast':
                ## Check for execution parameters required by BLAST method
                if e_values is None:
                    raise WorkflowError("You must specify at least one "
                                        "E-value.")
                ## Generate command for BLAST
                commands = _generate_blast_commands(output_dataset_dir,
                                                    input_fasta_fp,
                                                    reference_seqs_fp,
                                                    id_to_taxonomy_fp,
                                                    clean_otu_table_fp,
                                                    e_values)
                        
            ## Method is Mothur
            elif method == 'mothur':
                ## Check for execution parameters required by Mothur method
                if confidences is None:
                    raise WorkflowError("You must specify at least one "
                                        "confidence level.")
                ## Generate command for mothur
                commands = _generate_mothur_commands(output_dataset_dir,
                                                     input_fasta_fp,
                                                     reference_seqs_fp,
                                                     id_to_taxonomy_fp,
                                                     clean_otu_table_fp,
                                                     confidences)

            ## Method is RTAX
            elif method == 'rtax':
                ## Check for execution parameters required by RTAX method
                if rtax_modes is None:
                    raise WorkflowError("You must specify at least one mode "
                                        "to run RTAX in.")
                for mode in rtax_modes:
                    if mode not in ['single', 'paired']:
                        raise WorkflowError("Invalid rtax mode '%s'. Must be "
                                            "'single' or 'paired'." % mode)

                ## Generate command for rtax
                commands = _generate_rtax_commands(output_dataset_dir,
                                                   input_fasta_fp,
                                                   reference_seqs_fp,
                                                   id_to_taxonomy_fp,
                                                   clean_otu_table_fp,
                                                   rtax_modes,
                                                   read_1_seqs_fp,
                                                   read_2_seqs_fp)

            ## Unsupported method
            else:
                raise WorkflowError("Unrecognized or unsupported taxonomy "
                                    "assignment method '%s'." % method)

            # send command for current method to command handler
            for command in commands:
                start = time()

                # call_commands_serially needs a list of commands so here's a
                # length one commmand list.
                command_handler([command], status_update_callback, logger,
                                close_logger_on_success=False)
                end = time()
                logger.write('Time (s): %d\n\n' % (end - start))

    logger.close()

def _directory_check(output_dir, base_str, param_str):
    """ Checks to see if directories already exist from a previous run. """
    ## Save final and working output directory names
    final_dir = join(output_dir, base_str + param_str)
    working_dir = final_dir + '.tmp'
    ## Check if temp directory already exists (and delete if necessary)
    if isdir(working_dir):
        try:
            rmtree(working_dir)
        except OSError:
            raise WorkflowError("Temporary output directory already exists "
                                "(from a previous run perhaps) and cannot be "
                                "removed.")
    return final_dir, working_dir

def _generate_rdp_commands(output_dir, input_fasta_fp, reference_seqs_fp,
                           id_to_taxonomy_fp, clean_otu_table_fp, confidences,
                           rdp_max_memory=None):
    """ Build command strings for RDP method. """
    result = []
    for confidence in confidences:
        run_id = 'RDP, confidence: %s' % str(confidence)
        ## Get final and working directory names
        final_dir, working_dir = \
             _directory_check(output_dir, 'rdp_', str(confidence))
        ## Check if final directory already exists (skip iteration if it does)
        if isdir(final_dir):
            continue
        assign_taxonomy_command = \
                'assign_taxonomy.py -i %s -o %s -c %s -m rdp -r %s -t %s' % (
                input_fasta_fp, working_dir, str(confidence),
                reference_seqs_fp, id_to_taxonomy_fp)
        if rdp_max_memory is not None:
            assign_taxonomy_command += ' --rdp_max_memory %s' % rdp_max_memory
        result.append([('Assigning taxonomy (%s)' % run_id,
                      assign_taxonomy_command)])
        result.extend(_generate_taxa_processing_commands(working_dir,
                      input_fasta_fp, clean_otu_table_fp, run_id))
        ## Rename output directory
        result.append([('Renaming output directory (%s)' % run_id,
                      'mv %s %s' % (working_dir, final_dir))])
    return result

def _generate_blast_commands(output_dir, input_fasta_fp, reference_seqs_fp,
                             id_to_taxonomy_fp, clean_otu_table_fp, e_values):
    """ Build command strings for BLAST method. """
    result = []
    for e_value in e_values:
        run_id = 'BLAST, E: %s' % str(e_value)
        ## Get final and working directory names
        final_dir, working_dir = \
            _directory_check(output_dir, 'blast_', str(e_value))
        ## Check if final directory already exists (skip iteration if it does)
        if isdir(final_dir):
            continue
        assign_taxonomy_command = \
               'assign_taxonomy.py -i %s -o %s -e %s -m blast -r %s -t %s' % (
               input_fasta_fp, working_dir, str(e_value), reference_seqs_fp,
               id_to_taxonomy_fp)
        result.append([('Assigning taxonomy (%s)' % run_id,
                      assign_taxonomy_command)])
        result.extend(_generate_taxa_processing_commands(working_dir,
                      input_fasta_fp, clean_otu_table_fp, run_id))
        ## Rename output directory
        result.append([('Renaming output directory (%s)' % run_id,
                      'mv %s %s' % (working_dir, final_dir))])
    return result

def _generate_mothur_commands(output_dir, input_fasta_fp, reference_seqs_fp,
                              id_to_taxonomy_fp, clean_otu_table_fp,
                              confidences):
    """ Build command strings for Mothur method. """
    result = []
    for confidence in confidences:
        run_id = 'Mothur, confidence: %s' % str(confidence)
        ## Get final and working directory names
        final_dir, working_dir = \
            _directory_check(output_dir, 'mothur_', str(confidence))
        # Check if final directory already exists (skip iteration if it does)
        if isdir(final_dir):
            continue
        assign_taxonomy_command = \
               'assign_taxonomy.py -i %s -o %s -c %s -m mothur -r %s -t %s' % (
               input_fasta_fp, working_dir, str(confidence), reference_seqs_fp,
               id_to_taxonomy_fp)
        result.append([('Assigning taxonomy (%s)' % run_id,
                      assign_taxonomy_command)])
        result.extend(_generate_taxa_processing_commands(working_dir,
                      input_fasta_fp, clean_otu_table_fp, run_id))
        ## Rename output directory
        result.append([('Renaming output directory (%s)' % run_id,
                      'mv %s %s' % (working_dir, final_dir))])
    return result

def _generate_rtax_commands(output_dir, input_fasta_fp, reference_seqs_fp,
                            id_to_taxonomy_fp, clean_otu_table_fp, rtax_modes,
                            read_1_seqs_fp, read_2_seqs_fp):
    """ Build command strings for RTAX method. """
    result = []
    for mode in rtax_modes:
        run_id = 'RTAX, mode: %s' % mode
        ## Get final and working directory names
        final_dir, working_dir = \
                _directory_check(output_dir, 'rtax_', mode)
        ## Check if final directory already exists (skip iteration if it does)
        if isdir(final_dir):
            continue

        assign_taxonomy_command = \
                'assign_taxonomy.py -i %s -o %s -m rtax -r %s -t %s '\
                '--read_1_seqs_fp %s' % (input_fasta_fp, working_dir,
                reference_seqs_fp, id_to_taxonomy_fp, read_1_seqs_fp)
        if mode == 'paired':
            assign_taxonomy_command += ' --read_2_seqs_fp %s' % read_2_seqs_fp

        result.append([('Assigning taxonomy (%s)' % run_id,
                      assign_taxonomy_command)])
        result.extend(_generate_taxa_processing_commands(working_dir,
                      input_fasta_fp, clean_otu_table_fp, run_id))
        ## Rename output directory
        result.append([('Renaming output directory (%s)' % run_id,
                      'mv %s %s' % (working_dir, final_dir))])

    return result

def _generate_taxa_processing_commands(assigned_taxonomy_dir, input_fasta_fp,
                                       clean_otu_table_fp, run_id):
    """ Build command strings for adding and summarizing taxa commands. These 
        are used with every method. """
    taxa_assignments_fp = join(assigned_taxonomy_dir,
            splitext(basename(input_fasta_fp))[0] + '_tax_assignments.txt')
    otu_table_w_taxa_fp = join(assigned_taxonomy_dir,
            add_filename_suffix(clean_otu_table_fp, '_w_taxa'))
    add_md_command = [('Adding metadata (%s)' % run_id,
                       'add_metadata.py -i %s -o %s '
                       '--observation_mapping_fp %s --sc_separated taxonomy '
                       '--observation_header OTUID,taxonomy' %
                       (clean_otu_table_fp, otu_table_w_taxa_fp,
                        taxa_assignments_fp))]
    summarize_taxa_command = [('Summarizing taxa (%s)' % run_id,
                               'summarize_taxa.py -i %s -o %s' %
                               (otu_table_w_taxa_fp, assigned_taxonomy_dir))]

    return add_md_command, summarize_taxa_command
