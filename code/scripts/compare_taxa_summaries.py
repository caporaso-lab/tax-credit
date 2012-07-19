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

from os.path import join
from cogent.util.misc import create_dir
from qiime.parse import parse_taxa_summary_table
from qiime.util import (parse_command_line_parameters,
                        get_options_lookup,
                        make_option)

from taxcompare.compare_taxa_summaries import (add_filename_suffix,
        comparison_modes, compare_taxa_summaries, correlation_types,
        parse_sample_id_map)
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
        help='the two input taxa summary filepaths, comma-separated'),
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
     make_option('-s','--sample_id_map_fp', type='existing_filepath',
        help='Map of sample IDs from the first taxa summary file to compare '
        'against the sample IDs from the second taxa summary file. Each line '
        'should contain two sample IDs separated by a tab. More than one '
        'sample ID from the first taxa summary file may map to the same '
        'sample ID in the second taxa summary file if desired. This option '
        'only applies if the comparison mode is "paired". If not provided, '
        'only sample IDs that exist in both taxa summary files will be '
        'compared [default: %default]', default=None)
]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if len(opts.taxa_summary_fps) != 2:
        option_parser.error("Exactly two taxa summary files are required.")

    # Create the output dir if it doesn't already exist.
    try:
        create_dir(opts.output_dir)
    except:
        option_parser.error("Could not create or access output directory "
                            "specified with the -o option.")

    sample_id_map = None
    if opts.sample_id_map:
        sample_id_map = parse_sample_id_map(opts.sample_id_map)

    results = compare_taxa_summaries(
            parse_taxa_summary_table(open(opts.taxa_summary_fps[0], 'U')),
            parse_taxa_summary_table(open(opts.taxa_summary_fps[1], 'U')),
            opts.comparison_mode, correlation_type=opts.correlation_type,
            sample_id_map=sample_id_map)

    # Write out the sorted and filled taxa summaries, basing their
    # filenames on the original input filenames.
    for taxa_summary_fp, filled_taxa_summary in zip(opts.taxa_summary_fps,
                                                    results[:2]):
        taxa_summary_filled_fp = add_filename_suffix(taxa_summary_fp,
                                                     '_sorted_and_filled')
        filled_taxa_summary_f = open(join(opts.output_dir,
                                          taxa_summary_filled_fp), 'w')
        filled_taxa_summary_f.write(format_taxa_summary(filled_taxa_summary))
        filled_taxa_summary_f.close()

    # Write the correlation vector.
    corr_vec = open(join(opts.output_dir, 'all_samples_vs_expected.txt'), 'w')
    comp_matrix_f.write(results[2])
    comp_matrix_f.close()


if __name__ == "__main__":
    main()
