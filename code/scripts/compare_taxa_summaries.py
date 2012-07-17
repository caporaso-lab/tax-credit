#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout", "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"
__status__ = "Development"

from qiime.parse import parse_taxa_summary_table
from qiime.util import (parse_command_line_parameters,
                        get_options_lookup,
                        make_option)

from taxcompare.compare_taxa_summaries import (compare_taxa_summaries,
        comparison_modes, correlation_types)

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = """
"""
script_info['script_usage'] = []
script_info['script_usage'].append(("", "", ""))
script_info['output_description']= """
"""

script_info['required_options'] = [
    make_option('-i', '--taxa_summary_fps', type='existing_filepaths',
        help='the input taxa summary filepaths, comma-separated'),
    options_lookup['output_dir'],
    make_option('-m', '--comparison_mode', type='choice',
        choices=comparison_modes, help='the type of comparison to '
        'perform. "paired" will compare each sample in the taxa summary '
        'files. "expected" will compare each sample in the first taxa summary '
        'file to an expected taxa summary (specified in the second taxa '
        'summary file). If "expected", the second taxa summary file must '
        'contain only a single sample that all other samples will be compared '
        'to')
]
script_info['optional_options'] = [
    make_option('-c', '--correlation_type', type='choice',
        choices=correlation_types, help='the type of correlation '
        'measure to use [default: %default]', default='pearson'),
    make_option('-n', '--normalized_count', type='int',
        help='multiply all values in the taxa summary by this value. '
        'Useful to, for example, normalize frequency table to median '
        'seqs/sample count [default: %default]', default=1)
]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if len(opts.taxa_summary_fps) != 2:
        option_parser.error("Exactly two taxa summary files are required.")

    results = compare_taxa_summaries(
            parse_taxa_summary_table(open(opts.taxa_summary_fps[0], 'U')),
            parse_taxa_summary_table(open(opts.taxa_summary_fps[1], 'U')),
            opts.comparison_mode,
            correlation_type=opts.correlation_type,
            normalized_count=opts.normalized_count)


if __name__ == "__main__":
    main()
