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
            ## method is rdp
            if method == 'rdp':
                if rdp_id_to_taxonomy_fp is None:
                    raise WorkflowError("You must provide an ID to taxonomy "
                                        "map (formatted for RDP) filepath.")
                if confidences is None:
                    raise WorkflowError("You must specify at least one "
                                        "confidence level.")
                commands = _generate_rdp_commands(output_dataset_dir,
                                                  input_fasta_fp,
                                                  reference_seqs_fp,
                                                  rdp_id_to_taxonomy_fp,
                                                  clean_otu_table_fp,
                                                  confidences)

            ## method is blast
            else if method == 'blast':
                # check for necessary execution parameters for blast method
                if blast_id_to_taxonomy_fp is None:
                    raise WorkflowError("You must provide an ID to taxonomy "
                                        "map (formatted for BLAST) filepath.")
                if reference_seqs_fp is None:
                    raise WorkflowError("You must provide a reference sequences "
                                        "filepath.")
                if e_values is None:
                    raise WorkflowError("You must provide a maximum E-value "
                                        "for assignmnet using using BLAST.")

                # generate assign_taxonomy.py command parameters
                commands = _generate_blast_commands(output_dataset_dir,
                                                    input_fasta_fp,
                                                    reference_seqs_fp,
                                                    blast_id_to_taxonomy_fp,
                                                    clean_otu_table_fp,
                                                    e_values)

            ## method is mothur
            else if method == 'mothur':
                # check for necessary execution parameters for mothur method
                if mothur_id_to_taxonomy _fp is None:
                    raise WorkflowError("You must provide an ID to taxonomy "
                                        "map (formatted for mothur) filepath.")
                if reference_seqs_fp is None:
                    raise WorkflowError("You must provide a reference sequences "
                                        "filepath.")

                # generate mothur classify.seqs() command parameters
                commands = _generate_mothur_commands(output_dataset_dir,
                                                     input_fasta_fp,
                                                     reference_seqs_fp,
                                                     m_method, m_search, m_ksize, m_taxonomy,
                                                     m_procs, m_mistmatch, m_gapopen, m_save,
                                                     m_match, m_gapextend, m_numwanted, m_cutoff,
                                                     m_probs, m_iters)

            else:
                raise WorkflowError("Unrecognized or unsupported taxonomy "
                        "assignment method '%s'." % method)

            # call command_handler to execute commands
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
                                                         input_fasta_fp,
                                                         clean_otu_table_fp,
                                                         run_id))
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

def _generate_blast_commands(output_dir, input_fasta_fp, reference_seqs_fp,
                             blast_id_to_taxonomy_fp, clean_otu_table_fp,
                             e_values):
    result = []
    for e_value in e_values:
        run_id = 'BLAST, %s E' % str(e_value)
        assigned_taxonomy_dir = join(output_dir, 'blast_' + str(e_value))
        assign_taxonomy_command = \
                'assign_taxonomy.py -i %s -o %s -t %s -m blast -r %s -e %s' % (
                input_fasta_fp, output_dir, blast_id_to_taxonomy_fp, reference_seqs_fp, str(e_value))
        result.append([('Assigning taxonomy (%s)' % run_id, assign_taxonomy_command)])
        result.extend(_generate_taxa_processing_commands(assigned_taxonomy_dir,
                                                         input_fasta_fp,
                                                         clean_otu_table_fp,
                                                         run_id))
    return result

def _generate_mothur_commands():
    result = []
    for iter in range(1, iters):
        run_id = 'mothur, iteration %s' % str(iter)
        classify_seqs_dir = join(output_dir, 'mothur_' + str(iter))
        classify_seqs_command = \
                'classify.seqs(fasta=%s, reference=%s, method=%s, search=%s, ksize=%s, taxonomy=%s, processors=%s)' % ()
        result.append([('Classifying sequences (%s)' % run_id, classify_seqs_command)])
        result.extend(_generate_taxa_processing_commands(classify_seqs_dir,
                                                         input_fasta_fp,
                                                         clean_otu_table_fp,
                                                         run_id))
    return result
