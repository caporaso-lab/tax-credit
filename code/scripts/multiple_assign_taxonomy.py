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

from qiime.util import (parse_command_line_parameters, get_options_lookup,
                        make_option)

from taxcompare.multiple_assign_taxonomy import assign_taxonomy_multiple_times

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "Assigns taxonomy with multiple taxonomy assigners"
script_info['script_description'] = ""

script_info['script_usage'] = []
script_info['script_usage'].append(("", "", "%prog -h"))

script_info['output_description'] = ""

script_info['required_options'] = [
    make_option('-i', '--input_dir', type='existing_dirpath', help=''),
    make_option('-m', '--assignment_methods', type='string',
        help='Comma-separated list of taxon assignment methods to use, either '
        'blast, mothur, rdp, or rtax'),
    make_option('-r', '--reference_seqs_fp', type='existing_filepath',
        help='Path to reference sequences.  For assignment with blast, these '
        'are used to generate a blast database. For assignment with rdp, they '
        'are used as training sequences for the classifier')
]
script_info['optional_options'] = [
    make_option('--rdp_id_to_taxonomy_fp', type='existing_filepath',
        help='Path to tab-delimited file mapping sequences to assigned '
         'taxonomy. Each assigned taxonomy is provided as a '
         'semicolon-separated list. Each assigned taxonomy must be exactly 6 '
         'levels deep and each genus must be unique. REQUIRED when '
         'assignment_methods includes rdp [default: %default]', default=None),
    make_option('--blast_id_to_taxonomy_fp', type='existing_filepath',
        help='Path to tab-delimited file mapping sequences to assigned '
         'taxonomy. Each assigned taxonomy is provided as a '
         'semicolon-separated list. REQUIRED when assignment_methods includes '
         'blast [default: %default]', default=None),
    make_option('-c', '--confidences', type='string',
        help='Comma-separated list of minimum confidences to record an '
        'assignment, only used for rdp and mothur methods [default: %default]',
        default=None),
    make_option('-e', '--e_values', type='string',
        help='Comma-separated list of maximum e-values to record an '
        'assignment, only used for blast method [default: %default]',
        default=None),
    make_option('--rep_set_filename', type='string',
        help='[default: %default]', default='rep_set.fna'),
    make_option('--clean_otu_table_filename', type='string',
        help='[default: %default]', default='otu_table_mc2.biom'),
    make_option('-w', '--print_only', action='store_true',
        help='Print the commands but don\'t call them -- useful for debugging '
        '[default: %default]', default=False),
    make_option('-f', '--force', action='store_true',
        help='Force overwrite of existing taxonomy assignment directories '
        'under the input directory (note: only existing taxonomy assignment '
        'directories under input_dir will be removed). For example, if the '
        'directory "rdp_0.80" exists under input_dir/some_study/ and you run '
        'this script with rdp as one of the assignment methods and 0.80 as '
        'one of the confidence levels, this directory will be removed and '
        'recreated with the new output [default: %default]',
        default=False)
]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    assignment_methods = opts.assignment_methods.split(',')
    confidences = map(float, opts.confidences.split(','))
    e_values = map(float, opts.e_values.split(','))

    if opts.print_only:
        command_handler = print_commands
    else:
        command_handler = call_commands_serially

    if opts.verbose:
        status_update_callback = print_to_stdout
    else:
        status_update_callback = no_status_updates

    assign_taxonomy_multiple_times(opts.input_dir, assignment_methods,
        opts.reference_seqs_fp, opts.rep_set_filename,
        opts.clean_otu_table_filename,
        rdp_id_to_taxonomy_fp=opts.rdp_id_to_taxonomy_fp,
        blast_id_to_taxonomy_fp=opts.blast_id_to_taxonomy_fp,
        confidences=confidences, e_values=e_values,
        command_handler=command_handler,
        status_update_callback=status_update_callback, force=opts.force)


if __name__ == "__main__":
    main()
