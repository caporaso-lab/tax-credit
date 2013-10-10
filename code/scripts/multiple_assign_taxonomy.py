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
from qiime.workflow.util import (call_commands_serially, no_status_updates,
                                 print_commands, print_to_stdout)

from taxcompare.multiple_assign_taxonomy import (
        assign_taxonomy_multiple_times, split_input_str)

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "Assigns taxonomy with multiple taxonomy assigners"
script_info['script_description'] = """
This script provides a workflow to assign taxonomy using multiple taxonomy
assigners over multiple datasets.

It accepts a list of directories as input (one for each dataset), where each
directory contains the same set of files following a common naming convention.
At a minimum, each directory must have an OTU table without taxonomic
information and a fasta file containing the representative sequences that need
to be assigned taxonomy. Additional files may be needed depending on the
taxonomy assignment method; please refer to the script option documentation
below for more details.

The script creates (under the specified output directory) a subdirectory for
each input dataset directory. Under each of these directories, a subdirectory
is created for every possible combination of taxonomy assigners and parameters.
Each of these directories will contain an OTU table with taxonomic information,
the representative set's taxonomic assignments, and a set of taxa summary files
at varying taxonomic levels.

For example, if ran the script using the RDP classifier at 0.6 and 0.8
confidences over the S16S-1 input dataset, we'd get an output directory that
looks like:

output_dir
  S16S-1
    rdp_0.6
      otu_table.biom
      ...
    rdp_0.8
      otu_table.biom
      ...

The script supports efficient reruns, meaning that it will not reprocess a
taxonomy assigner/parameter combo if the expected output subdirectory already
exists. This allows the script to be rerun in an incremental fashion (e.g. if
the script fails or we decide to add a new assigner or parameter, we don't want
to have to rerun everything from scratch).
"""

script_info['script_usage'] = []
script_info['script_usage'].append(("Multiple RDP and mothur assignments",
"This example shows how to assign taxonomy using RDP and mothur, at 0.6 and "
"0.8 confidence levels, over the S16S-1 and S16S-2 datasets. This command "
"will result in a total of eight subdirectories (four for each dataset).",
"%prog -i S16S-1,S16S-2 -o example_1_output -m rdp,mothur -r gg_97_otus_4feb2011.fasta -t greengenes_tax.txt -c 0.6,0.8"))

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
    make_option('--uclust_min_consensus_fractions', type='string',
        help='Comma-separated list of floats indicating minimum consensus '
        'fractions. Each value specifies the minimum fraction of database '
        'hits that must have a specific taxonomic assignment to assign that '
        'taxonomy to a query, only used for uclust method [default: %default]',
        default=None),
    make_option('--uclust_similarities', type='string',
        help='Comma-separated list of floats indicating the minimum percent '
        'similarities to consider a database match a hit, only used for '
        'uclust method [default: %default]', default=None),
    make_option('--uclust_max_accepts', type='string',
        help='Comma-separated list of integers indicating the number of '
        'database hits to consider when making an assignment, only used for '
        'uclust method [default: %default]', default=None),
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
    make_option('--rtax_read_id_regexes', type='string',
        help='Comma-separated list of regular expressions that will be passed '
        'to RTAX for parsing read 1 IDs. The number (and order) of the '
        'supplied regexes must match the number (and order) of input datasets '
        'supplied via -i/--input_dirs [default: %default]', default=None),
    make_option('--rtax_amplicon_id_regexes', type='string',
        help='Comma-separated list of regular expressions that will be passed '
        'to RTAX for parsing amplicon IDs. The number (and order) of the '
        'supplied regexes must match the number (and order) of input datasets '
        'supplied via -i/--input_dirs [default: %default]', default=None),
    make_option('--rtax_header_id_regexes', type='string',
        help='Comma-separated list of regular expressions that will be passed '
        'to RTAX for parsing header IDs. The number (and order) of the '
        'supplied regexes must match the number (and order) of input datasets '
        'supplied via -i/--input_dirs [default: %default]', default=None),
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

    confidences = split_input_str(opts.confidences)
    e_values = split_input_str(opts.e_values)
    rtax_modes = split_input_str(opts.rtax_modes, map_fn=str)
    uclust_min_consensus_fractions = \
            split_input_str(opts.uclust_min_consensus_fractions)
    uclust_similarities = split_input_str(opts.uclust_similarities)
    uclust_max_accepts = split_input_str(opts.uclust_max_accepts, map_fn=int)

    rtax_read_id_regexes = split_input_str(opts.rtax_read_id_regexes,
                                           map_fn=str)
    rtax_amplicon_id_regexes = split_input_str(opts.rtax_amplicon_id_regexes,
                                               map_fn=str)
    rtax_header_id_regexes = split_input_str(opts.rtax_header_id_regexes,
                                             map_fn=str)

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
        e_values=e_values, rtax_modes=rtax_modes,
        uclust_min_consensus_fractions=uclust_min_consensus_fractions,
        uclust_similarities=uclust_similarities,
        uclust_max_accepts=uclust_max_accepts,
        input_fasta_filename=opts.input_fasta_filename,
        clean_otu_table_filename=opts.clean_otu_table_filename,
        read_1_seqs_filename=opts.read_1_seqs_filename,
        read_2_seqs_filename=opts.read_2_seqs_filename,
        rtax_read_id_regexes=rtax_read_id_regexes,
        rtax_amplicon_id_regexes=rtax_amplicon_id_regexes,
        rtax_header_id_regexes=rtax_header_id_regexes,
        rdp_max_memory=opts.rdp_max_memory,
        command_handler=command_handler,
        status_update_callback=status_update_callback, force=opts.force)


if __name__ == "__main__":
    main()
