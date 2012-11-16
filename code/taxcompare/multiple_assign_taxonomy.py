#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"
__status__ = "Development"

"""Contains functions used in the multiple_assign_taxonomy.py script."""

from os import makedirs, rename
from os.path import basename, isdir, join, normpath, split, splitext
from shutil import rmtree
from qiime.util import add_filename_suffix
from qiime.workflow import (call_commands_serially, generate_log_fp,
                            no_status_updates, print_commands, print_to_stdout,
                            WorkflowError, WorkflowLogger)

def assign_taxonomy_multiple_times(input_dirs, output_dir, assignment_methods,
        reference_seqs_fp, input_fasta_filename, clean_otu_table_filename,
        id_to_taxonomy_fp=None, confidences=None, e_values=None,
        command_handler=call_commands_serially, rdp_max_memory=None,
        status_update_callback=print_to_stdout, force=False,
        read_1_seqs_fp=None, read_2_seqs_fp=None):
    ## Check if temp output directory exists
    try:
        makedirs(output_dir)
    except OSError:
        if not force:
            raise WorkflowError("Output directory '%s' already exists. Please "
                    "choose a different directory, or force overwrite with -f."
                    % output_dir)

    # Check for inputs that are universally required
    if assignment_methods is None:
        raise WorkflowError("You must specify at least one method:" 
                            "'rdp', 'blast', 'mothur', or 'rtax'.")
    if input_fasta_filename is None:
        raise WorkflowError("You must provide an input fasta filename.")
    if clean_otu_table_filename is None:
        raise WorkflowError("You must provide a clean otu table filename.")
    if id_to_taxonomy_fp is None:
        raise WorkflowError("You must provide an ID to taxonomy map filename.")
    
    logger = WorkflowLogger(generate_log_fp(output_dir))

    for input_dir in input_dirs:
        # Make sure the input dataset directory exists.
        if not isdir(input_dir):
            raise WorkflowError("The input directory '%s' does not exist." %
                                input_dir)

        input_dir_name = split(normpath(input_dir))[1]
        output_dataset_dir = join(output_dir, input_dir_name)
        input_fasta_fp = join(input_dir, input_fasta_filename)
        clean_otu_table_fp = join(input_dir, clean_otu_table_filename)

        logger.write("\nCreating output subdirectory '%s' if it doesn't "
                     "already exist.\n" % output_dataset_dir)
        try:
            makedirs(output_dataset_dir)
        except OSError:
            # It already exists, which is okay since we already know we are in
            # 'force' mode from above.
            pass

        for method in assignment_methods:
            # method is RDP
            if method == 'rdp':
                # check for execution parameters required by RDP method
                if confidences is None:
                    raise WorkflowError("You must specify at least one "
                                        "confidence level.")
                # generate command for RDP
                commands = _generate_rdp_commands(output_dataset_dir,
                                                  input_fasta_fp,
                                                  reference_seqs_fp,
                                                  id_to_taxonomy_fp,
                                                  clean_otu_table_fp,
                                                  confidences,
                                                  rdp_max_memory=rdp_max_memory)
                        
            # method is BLAST
            elif method == 'blast':
                # check for execution parameters required by BLAST method
                if e_values is None:
                    raise WorkflowError("You must specify at least one "
                                        "E value.")
                # generate command for BLAST
                commands = _generate_blast_commands(output_dataset_dir,
                                                    input_fasta_fp,
                                                    reference_seqs_fp,
                                                    id_to_taxonomy_fp,
                                                    clean_otu_table_fp,
                                                    e_values)
                        
            # method is Mothur
            elif method == 'mothur':
                # check for execution parameters required by Mothur method
                if confidences is None:
                    raise WorkflowError("You must specify at least one "
                                        "confidence level.")
                # generate command for mothur
                commands = _generate_mothur_commands(output_dataset_dir,
                                                     input_fasta_fp,
                                                     reference_seqs_fp,
                                                     id_to_taxonomy_fp,
                                                     clean_otu_table_fp,
                                                     confidences)

            # method is RTAX
            elif method == 'rtax':
                # check for execution parameters required by Mothur method
                if read_1_seqs_fp is None:
                    raise WorkflowError("You must specify a file containing the first "
                                        "read from pair-end sequencing.")
                # generate command for mothur
                commands = _generate_rtax_commands(output_dataset_dir,
                                                   input_fasta_fp,
                                                   reference_seqs_fp,
                                                   id_to_taxonomy_fp,
                                                   clean_otu_table_fp,
                                                   read_1_seqs_fp,
                                                   read_2_seqs_fp=read_2_seqs_fp)

            # unsupported method
            else:
                raise WorkflowError("Unrecognized or unsupported taxonomy "
                        "assignment method '%s'." % method)
            # send command for current method to command handler
            command_handler(commands, status_update_callback, logger,
                            close_logger_on_success=False)
    logger.close()

def _generate_rdp_commands(output_dir, input_fasta_fp, reference_seqs_fp,
                           id_to_taxonomy_fp, clean_otu_table_fp, confidences,
                           rdp_max_memory=None):
    result = []
    for confidence in confidences:
        run_id = 'RDP, %s confidence' % str(confidence)
        # Save final and working output directory names
        final_dir = join(output_dir, 'rdp_' + str(confidence))
        working_dir = final_dir + '.tmp'
        # Check if temp directory already exists (and delete if necessary)
        if isdir(working_dir):
            try:
                rmtree(working_dir)
            except OSError:
                raise WorkflowError("Temporary output directory exists (from a "
                                    "previous run perhaps) and cannot be removed.")
        # Check if final directory already exists (skip iteration if it does)
        if isdir(final_dir):
            continue
        assign_taxonomy_command = \
                'assign_taxonomy.py -i %s -o %s -c %s -m rdp -r %s -t %s' % (
                input_fasta_fp, working_dir, str(confidence), reference_seqs_fp,
                id_to_taxonomy_fp)
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
    result = []
    for e in e_values:
        run_id = 'BLAST, E %s' % str(e)
        # Save final and working output directory names
        final_dir = join(output_dir, 'blast_' + str(e))
        working_dir = final_dir + '.tmp'
        # Check if temp directory already exists (and delete if necessary)
        if isdir(working_dir):
            try:
                rmtree(working_dir)
            except OSError:
                raise WorkflowError("Temporary output directory exists (from a "
                                    "previous run perhaps) and cannot be removed.")
        # Check if final directory already exists (skip iteration if it does)
        if isdir(final_dir):
            continue
        assign_taxonomy_command = \
                'assign_taxonomy.py -i %s -o %s -e %s -m blast -r %s -t %s' % (
                input_fasta_fp, working_dir, str(e), reference_seqs_fp,
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
                              id_to_taxonomy_fp, clean_otu_table_fp, confidences):
    result = []
    for confidence in confidences:
        run_id = 'Mothur, %s confidence' % str(confidence)
        # Save final and working output directory names
        final_dir = join(output_dir, 'mothur_' + str(confidence))
        working_dir = final_dir + '.tmp'
        # Check if temp directory already exists (and delete if necessary)
        if isdir(working_dir):
            try:
                rmtree(working_dir)
            except OSError:
                raise WorkflowError("Temporary output directory exists (from a "
                                    "previous run perhaps) and cannot be removed.")
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
                            id_to_taxonomy_fp, clean_otu_table_fp,
                            read_1_seqs_fp, read_2_seqs_fp=None):
    result = []
    run_id = 'RTAX'
    for i in range(1,2):
        ## For single-end reads
        if i is 1:
            run_id = 'RTAX, single-end'
            final_dir = join(output_dir, 'rtax_single')
        ## For paired-end reads
        else:
            run_id = 'RTAX, paired-end'
            final_dir = join(output_dir, 'rtax_paired')
        # Save working output directory name
        working_dir = final_dir + '.tmp'
        # Check if temp directory already exists (and delete if necessary)
        if isdir(working_dir):
            try:
                rmtree(working_dir)
            except OSError:
                raise WorkflowError("Temporary output directory exists (from a "
                                    "previous run perhaps) and cannot be removed.")
        # Check if final directory already exists (skip iteration if it does)
        if isdir(final_dir):
            continue
        assign_taxonomy_command = \
                'assign_taxonomy.py -i %s -o %s -m rtax -r %s -t %s '\
                '--read_1_seqs_fp %s' % (
                input_fasta_fp, working_dir, reference_seqs_fp, id_to_taxonomy_fp,
                read_1_seqs_fp)
        ## Append second read parameter for paired end
        if i is 2:
            assigned_taxonomy_command += ' --read_2_seqs_fp %s' % read_2_seqs_fp
        result.append([('Assigning taxonomy (%s)' % run_id,
                      assign_taxonomy_command)])
        result.extend(_generate_taxa_processing_commands(working_dir,
                      input_fasta_fp, clean_otu_table_fp, run_id))
        ## Rename output directory
        result.append([('Renaming output directory (%s)' % run_id,
                      'mv %s %s' % (working_dir, final_dir))])
        ## Break execution if no second read parameter is provided
        if read_2_seqs_fp is None:
            return result

    return result

def _generate_taxa_processing_commands(assigned_taxonomy_dir, input_fasta_fp,
                                       clean_otu_table_fp, run_id):
    taxa_assignments_fp = join(assigned_taxonomy_dir,
            splitext(basename(input_fasta_fp))[0] + '_tax_assignments.txt')
    otu_table_w_taxa_fp = join(assigned_taxonomy_dir,
            add_filename_suffix(clean_otu_table_fp, '_w_taxa'))

    add_taxa_command = [('Adding taxa (%s)' % run_id,
                        'add_taxa.py -i %s -o %s -t %s' %
                        (clean_otu_table_fp, otu_table_w_taxa_fp,
                         taxa_assignments_fp))]
    summarize_taxa_command = [('Summarizing taxa (%s)' % run_id,
                              'summarize_taxa.py -i %s -o %s' %
                              (otu_table_w_taxa_fp, assigned_taxonomy_dir))]

    return add_taxa_command, summarize_taxa_command
