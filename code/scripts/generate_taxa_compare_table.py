#!/usr/bin/env python

__author__ = "Kyle Patnode"
__copyright__ = "Copyright 2012, The QIIME Project"
__credits__ = ["Kyle Patnode", "Jai Ram Rideout", "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Kyle Patnode"
__email__ = "kpatnode1@gmail.com"
__status__ = "Development"

from itertools import izip
from qiime.util import parse_command_line_parameters, get_options_lookup, make_option, create_dir
from cogent.util.misc import create_dir
from taxcompare.generate_taxa_compare_table import generate_taxa_compare_table

options_lookup = get_options_lookup()

script_info={}
script_info['brief_description']="""Walks a file tree from a given root and compares the taxa summaries found."""
script_info['script_description'] = ""

script_info['script_usage']=[]

script_info['script_usage'].append(("", "", "%prog -h"))

script_info['output_description']="""A tab-delimited table showing Pearson's and Spearman's correalation between the expected (in the key files) and the actual (found within the root). There is a file for every level compared."""
script_info['required_options']=[\
 make_option('-r', '--root_dir',type="existing_dirpath",
        help='Path to the root of the output from multiple_assign_taxonomy.py'),\

 make_option('-k', '--key_dir',type="existing_dirpath",
        help='Path to file containing the expected results of the files in the output directory.')]

script_info['optional_options']=[\
 make_option('--output_dir', type="new_dirpath",
        help='Prefix for output files. If None is given, will print output.',
        default = None),\

 make_option('--levels', type="string",
        help='Comma-separated list of multiple_assign_taxonomy output levels to analyze. Numbers between 2 and 6 inclusive.'+\
        '[default: %default]',
        default = '2,3,4,5,6'),\

 make_option('--force', action='store_true',
	help='If true forces output directory creation even if the directory already exists.'+\
	'[default: %default]',
	default = False)]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if(opts.output_dir):
        create_dir(opts.output_dir, fail_on_exist=opts.force)

    levels = map(int, opts.levels.split(','))

    results = generate_taxa_compare_table(opts.root_dir, opts.key_dir, levels)

    for level, out_file in izip(levels, results):
        if(opts.output_dir):
            with open(opts.output_dir + '_L' + str(level) + '.txt', 'w') as f:
                datasets = sorted(out_file.keys())
                methods = set()
                for m in out_file.itervalues():
                    methods|= set(m.keys())
                methods = sorted(list(methods))
                f.write('P_S\t')
                f.write('\t'.join(methods)+'\n')
                for dataset in datasets:
                    f.write(dataset+'\t')
                    for method in methods:
                        try:
                            f.write(out_file[dataset][method][0] + "_" + out_file[dataset][method][1]+'\t')
                        except KeyError:
                            f.write('N/A'+'\t')#Don't have data for that set/method
                    f.write('\n')
        else:
            #Terminal output code here if wanted
            pass

if __name__ == '__main__':
    main()
