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

from qiime.util import (parse_command_line_parameters, get_options_lookup,
                        make_option)
from qiime.workflow import call_commands_serially, no_status_updates

from taxcompare.multiple_assign_taxonomy import assign_taxonomy_multiple_times

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "Assigns taxonomy with multiple taxonomy assigners"
script_info['script_description'] = """
This script provides a workflow to assign taxonomy using multiple taxonomy
assigners over multiple datasets.

It accepts a list of directories as input (one for each dataset), where each
directory contains the same set of files following a common naming convention.
At a minimum, each directory must have an 
"""

script_info['script_usage'] = []
script_info['script_usage'].append(("", "", "%prog -h"))

script_info['output_description'] = ""

script_info['required_options'] = [
    make_option('-i', '--input_dirs', type='string',
        help='Comma-separated list of input dataset directories'),
    options_lookup['output_dir'],
    make_option('-m', '--assignment_methods', type='string',
        help='Comma-separated list of taxon assignment methods to use, either '
        'blast, mothur, rdp, or rtax'),
    make_option('-r', '--reference_seqs_fp', type='existing_filepath',
        help='Path to reference sequences.  For assignment with blast, these '
        'are used to generate a blast database. For assignment with rdp and '
        'mothur, they are used as training sequences for the classifier'),
    make_option('-t', '--id_to_taxonomy_fp', type='existing_filepath',
        help='Path to tab-delimited file mapping sequences to assigned '
        'taxonomy. Each assigned taxonomy is provided as a '
        'semicolon-separated list.')
]
script_info['optional_options'] = [
    make_option('-c', '--confidences', type='string',
        help='Comma-separated list of minimum confidences to record an '
        'assignment, only used for rdp and mothur methods [default: %default]',
        default=None),
    make_option('-e', '--e_values', type='string',
        help='Comma-separated list of maximum e-values to record an '
        'assignment, only used for blast method [default: %default]',
        default=None),
    make_option('-x', '--rtax_modes', type='string',
        help='Comma-separated list of modes to run rtax in. Modes can be '
        'either "single" or "paired". If paired, paired-end reads must be '
        'available (i.e. forward and reverse reads, demultiplexed) '
        '[default: %default]', default=None),
    make_option('--input_fasta_filename', type='string',
        help='Name of fasta file containing sequences to receive taxonomy '
        'assignment. Must exist under each input dataset directory '
        '[default: %default]', default='rep_set.fna'),
    make_option('--clean_otu_table_filename', type='string',
        help='Name of OTU table BIOM file that will have taxonomic '
        'information added to it. Must exist under each input dataset '
        'directory [default: %default]', default='otu_table_mc2.biom'),
    make_option('--read_1_seqs_filename', type='string',
        help='Name of fasta file containing the first read from paired-end '
        'sequencing, prior to OTU clustering (used for RTAX only). Must exist '
        'under each input dataset directory if rtax is specified as an '
        'assignment method [default: %default]', default='seqs1.fna'),
    make_option('--read_2_seqs_filename', type='string',
        help='Name of fasta file containing a second read from paired-end '
        'sequencing, prior to OTU clustering (used for RTAX only). Must exist '
        'under each input dataset directory if rtax is specified as an '
        'assignment method and the "paired" rtax mode is specified '
        '[default: %default]', default='seqs2.fna'),
    make_option('--rdp_max_memory', type='string',
        help='Maximum memory allocation, in MB, for JVM when using the rdp '
        'method. Increase for large training sets [default: %default]',
        default=4000),
    make_option('-w', '--print_only', action='store_true',
        help='Print the commands but don\'t call them -- useful for debugging '
        '[default: %default]', default=False),
    make_option('-f', '--force', action='store_true',
        help='Force overwrite of existing output directory (note: existing '
        'files in output_dir will not be removed). Subdirectories that '
        'already exist will be skipped (they are assumed to have been '
        'successfully created during a previous run of the script) '
        '[default: %default]', default=False)
]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    input_dirs = opts.input_dirs.split(',')
    assignment_methods = opts.assignment_methods.split(',')

    confidences = opts.confidences
    if confidences is not None:
        confidences = map(float, opts.confidences.split(','))

    e_values = opts.e_values
    if e_values is not None:
        e_values = map(float, opts.e_values.split(','))

    rtax_modes = opts.rtax_modes
    if rtax_modes is not None:
        rtax_modes = opts.rtax_modes.split(',')

    if opts.print_only:
        command_handler = print_commands
    else:
        command_handler = call_commands_serially

    if opts.verbose:
        status_update_callback = print_to_stdout
    else:
        status_update_callback = no_status_updates

    assign_taxonomy_multiple_times(input_dirs, opts.output_dir,
        assignment_methods, opts.reference_seqs_fp,
        opts.id_to_taxonomy_fp, confidences=confidences,
        e_values=e_values, rtax_modes=opts.rtax_modes,
        input_fasta_filename=opts.input_fasta_filename,
        clean_otu_table_filename=opts.clean_otu_table_filename,
        read_1_seqs_filename=opts.read_1_seqs_filename,
        read_2_seqs_filename=opts.read_2_seqs_filename,
        rdp_max_memory=opts.rdp_max_memory,
        command_handler=command_handler,
        status_update_callback=status_update_callback, force=opts.force)


if __name__ == "__main__":
    main()
