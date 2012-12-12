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

from time import time
from os import makedirs
from os.path import basename, isdir, join, normpath, split, splitext
from qiime.util import add_filename_suffix
from qiime.workflow import (call_commands_serially, generate_log_fp,
                            no_status_updates, print_commands, print_to_stdout,
                            WorkflowError, WorkflowLogger)

def assign_taxonomy_multiple_times(input_dirs, output_dir, assignment_methods,
        reference_seqs_fp, input_fasta_filename, clean_otu_table_filename,
        id_to_taxonomy_fp=None, confidences=None, e_values=None,
        command_handler=call_commands_serially,
        status_update_callback=print_to_stdout, force=False):
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
                            "'rdp', 'blast', or 'mothur'.")
    if input_fasta_filename is None:
        raise WorkflowError("You must provide an input fasta filename.")
    if clean_otu_table_filename is None:
        raise WorkflowError("You must provide a clean otu table filename.")
    if id_to_taxonomy_fp is None:
        raise WorkflowError("You must provide an ID to taxonomy map filename.")
    
    logger = WorkflowLogger(generate_log_fp(output_dir))
    time_results=[]

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
                                                  confidences)
                        
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
                                                     reference_seqs,
                                                     id_to_taxonomy_fp,
                                                     clean_otu_table_fp,
                                                     confidences)
            # unsupported method
            else:
                raise WorkflowError("Unrecognized or unsupported taxonomy "
                        "assignment method '%s'." % method)

            # send command for current method to command handler
            for command in commands:
                #call_commands_serially needs a list of commands so here's a length one commmand list.
                c = list()
                c.append(command)
                start = time()
                command_handler(c, status_update_callback, logger,
                                close_logger_on_success=False)
                end = time()
                input_file = command[0][1].split()[command[0][1].split().index('-i')+1].split('/')[-2]
                if 'Assigning' in command[0][0]:
                    time_results.append((input_file, ' '.join(command[0][0].split()[2:]), end-start))

    # removes and writes out the title we initialized with earlier
    logger.write('\n\nAssignment times (seconds):\n')
    for t in time_results:
        # write out each time result as (method, params)\ttime (seconds)
        #First clean up the output
        method, param = t[1].split(', ')
        method = method.lstrip('(')
        param = param.rstrip(')')

        logger.write('%s\t%s\t%s\t%s\n' % (t[0], method, param, str(t[2])))

    logger.close()

def _generate_rdp_commands(output_dir, input_fasta_fp, reference_seqs_fp,
                           id_to_taxonomy_fp, clean_otu_table_fp, confidences):
    result = []
    for confidence in confidences:
        run_id = 'RDP, %s confidence' % str(confidence)
        assigned_taxonomy_dir = join(output_dir, 'rdp_' + str(confidence))
        assign_taxonomy_command = \
                'assign_taxonomy.py -i %s -o %s -c %s -m rdp -r %s -t %s' % (
                input_fasta_fp, assigned_taxonomy_dir, str(confidence),
                reference_seqs_fp, id_to_taxonomy_fp)
        result.append([('Assigning taxonomy (%s)' % run_id,
                       assign_taxonomy_command)])
        result.extend(_generate_taxa_processing_commands(assigned_taxonomy_dir,
            input_fasta_fp, clean_otu_table_fp, run_id))
    return result

def _generate_blast_commands(output_dir, input_fasta_fp, reference_seqs_fp,
                             id_to_taxonomy_fp, clean_otu_table_fp, e_values):
    result = []
    for e in e_values:
        run_id = 'BLAST, E %s' % str(e)
        assigned_taxonomy_dir = join(output_dir, 'blast_' + str(e))
        assign_taxonomy_command = \
                'assign_taxonomy.py -i %s -o %s -e %s -m blast -r %s -t %s' % (
                input_fasta_fp, assigned_taxonomy_dir, str(e),
                reference_seqs_fp, id_to_taxonomy_fp)
        result.append([('Assigning taxonomy (%s)' % run_id,
                       assign_taxonomy_command)])
        result.extend(_generate_taxa_processing_commands(assigned_taxonomy_dir,
            input_fasta_fp, clean_otu_table_fp, run_id))
    return result

def _generate_mothur_commands(output_dir, input_fasta_fp, reference_seqs_fp,
                              id_to_taxonomy_fp, clean_otu_table_fp, confidences):
    result = []
    for confidence in confidences:
        run_id = 'Mothur, %s confidence' % str(confidence)
        assigned_taxonomy_dir = join(output_dir, 'mothur_%s' % str(confidence))
        assign_taxonomy_command = \
                'assign_taxonomy.py -i %s -o %s -c %s -m mothur -r %s -t %s' % (
                input_fasta_fp, assigned_taxonomy_dir, str(confidence), 
                reference_seqs_fp, id_to_taxonomy_fp)
        result.append([('Assigning taxonomy (%s)' % run_id,
                       assign_taxonomy_command)])
        result.extend(_generate_taxa_processing_commands(assigned_taxonomy_dir,
            input_fasta_fp, clean_otu_table_fp, run_id))
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
