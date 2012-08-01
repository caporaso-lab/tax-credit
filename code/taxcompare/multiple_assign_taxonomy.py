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

from os import makedirs
from os.path import basename, isdir, join, normpath, split, splitext
from qiime.util import add_filename_suffix
from qiime.workflow import (call_commands_serially, generate_log_fp,
                            no_status_updates, print_commands, print_to_stdout,
                            WorkflowError, WorkflowLogger)

def assign_taxonomy_multiple_times(input_dirs, output_dir, assignment_methods,
        reference_seqs_fp, input_fasta_filename, clean_otu_table_filename,
        rdp_id_to_taxonomy_fp=None, blast_id_to_taxonomy_fp=None,
        confidences=None, e_values=None,
        command_handler=call_commands_serially,
        status_update_callback=print_to_stdout, force=False):
    try:
        makedirs(output_dir)
    except OSError:
        if not force:
            raise WorkflowError("Output directory '%s' already exists. Please "
                    "choose a different directory, or force overwrite with -f."
                    % output_dir)

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
            if method == 'rdp':
                if rdp_id_to_taxonomy_fp is None:
                    raise WorkflowError("You must provide an ID to taxonomy "
                                        "map (formatted for RDP) filepath.")
                if confidences is None:
                    raise WorkflowError("You must specify at least one "
                                        "confidence level.")
                commands = _generate_rdp_commands(output_dataset_dir,
                        input_fasta_fp, reference_seqs_fp,
                        rdp_id_to_taxonomy_fp, clean_otu_table_fp, confidences)
            else:
                raise WorkflowError("Unrecognized or unsupported taxonomy "
                        "assignment method '%s'." % method)

            command_handler(commands, status_update_callback, logger,
                            close_logger_on_success=False)
    logger.close()

def _generate_rdp_commands(output_dir, input_fasta_fp, reference_seqs_fp,
                           rdp_id_to_taxonomy_fp, clean_otu_table_fp,
                           confidences):
    result = []
    for confidence in confidences:
        run_id = 'RDP, %s confidence' % str(confidence)
        assigned_taxonomy_dir = join(output_dir, 'rdp_' + str(confidence))
        assign_taxonomy_command = \
                'assign_taxonomy.py -i %s -o %s -c %s -m rdp -r %s -t %s' % (
                input_fasta_fp, assigned_taxonomy_dir, str(confidence),
                reference_seqs_fp, rdp_id_to_taxonomy_fp)
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
