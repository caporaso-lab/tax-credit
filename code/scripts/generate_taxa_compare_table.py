#!/usr/bin/env python

__author__ = "Kyle Patnode"
__copyright__ = "Copyright 2012, The QIIME Project"
__credits__ = ["Kyle Patnode", "Jai Ram Rideout", "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Kyle Patnode"
__email__ = "kpatnode1@gmail.com"
__status__ = "Development"

from os.path import join
from itertools import izip
from qiime.util import parse_command_line_parameters, get_options_lookup, make_option, create_dir
from taxcompare.generate_taxa_compare_table import generate_taxa_compare_table, format_output

options_lookup = get_options_lookup()

script_info={}
script_info['brief_description']="""Walks a file tree from a given root and compares the taxa summaries found to the given keys."""
script_info['script_description'] = """Contains code to search a file tree for OTU files created by multiple_assign_taxonomy.py
and compare them to the OTU files containing expected compositions found in the key_dir. Comparison 
uses compare_taxa_summaries and gets both Pearson and Spearman correlations using a paired analysis."""

script_info['script_usage']=[]

script_info['script_usage'].append(("Sample Usage with user-defined levels:", "When given levels by the user, "
"generate_taxa_compare_table.py will only analyze those levels of output from multiple_assign_taxonomy.py. "
"This command will only compare levels 2, 4, and 5: ",
"%prog -r root_of_directory -k directory_containing_only_key_files -o output_dir -l 2,4,5"))

script_info['output_description']="""A tab-delimited table showing Pearson's and Spearman's correalation between the expected (in the key files) and the actual (found within the root). There is a file for every level compared."""
script_info['required_options']=[
 make_option('-r', '--root_dir',type="existing_dirpath",
        help='Path to the root of the output from multiple_assign_taxonomy.py'),

 make_option('-k', '--key_dir',type="existing_dirpath",
        help='Path to file containing the expected results of the files in the output directory.'),

 make_option('-o', '--output_dir', type="new_dirpath",
        help='Directory to output files.')]

script_info['optional_options']=[
 make_option('-l', '--levels', type="string",
        help='Comma-separated list of multiple_assign_taxonomy output levels to analyze. Numbers between 2 and 6 inclusive.'
        '[default: %default]',
        default = '2,3,4,5,6'),

 make_option('-s', '--separator', type="string",
        help='Sets the string separator used to split up Pearson and Spearman coefficients in the output.'
        '[default: %default]',
        default = ',')]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    create_dir(opts.output_dir, fail_on_exist=False)

    levels = map(int, opts.levels.split(','))

    results = generate_taxa_compare_table(opts.root_dir, opts.key_dir, levels)
    results = format_output(results, opts.separator)

    for level in levels:
        with open(join(opts.output_dir, 'compare_table_L' + str(level) + '.txt'), 'w') as f:
            f.writelines(results[level])

if __name__ == '__main__':
    main()
