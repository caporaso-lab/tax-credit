#!/usr/bin/env python

__author__ = "Kyle Patnode"
__copyright__ = "Copyright 2012, The QIIME Project"
__credits__ = ["Kyle Patnode", "Jai Ram Rideout", "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Kyle Patnode"
__email__ = "kpatnode1@gmail.com"
__status__ = "Development"

from os import walk
from qiime.parse import parse_taxa_summary_table
from qiime.compare_taxa_summaries import compare_taxa_summaries

assignment_method_choices = ['rdp','blast','rtax','mothur','tax2tree']

def get_key_files(directory):
    key_fps = {}
    for key_file in walk(directory).next()[2]:
        experiment = key_file.split('_')[0].capitalize()
        key_fps[experiment] = directory + key_file
    return key_fps

def get_coefficients(run_fp, key_fp):
    run = parse_taxa_summary_table(open(run_fp, 'U'))
    key = parse_taxa_summary_table(open(key_fp, 'U'))

    pearson_compare = compare_taxa_summaries(run, key, 'paired', 'pearson')
    spearman_compare = compare_taxa_summaries(run, key, 'paired', 'spearman')

    pearson_coeff = pearson_compare[2].split('\n')[-2].split()[0]
    spearman_coeff = spearman_compare[2].split('\n')[-2].split()[0]

    return pearson_coeff, spearman_coeff

def generate_taxa_compare_table(root, key_directory, levels):
    key_fps = get_key_files(key_directory)

    results = []
    for i in range(len(levels)):
        results.append(list())

    for(path, dirs, files) in walk(root, topdown=True):
        dirs.sort()
        files.sort()
        for choice in assignment_method_choices:
            if choice in path:
                experiment = path.split('/')[-2].rstrip('-123').capitalize()
                for f in files:
                    if 'otu_table_mc2_w_taxa_L' in f:
                        name = path.split('/')[-2].capitalize()
                        level = int(f[-5])-2
                        try:
                            pearson_coeff, spearman_coeff = get_coefficients(path+'/'+f, key_fps[experiment])
                        except ValueError:#compare_taxa_summaries couldn't find a match between the 2
                        #Likely due to mismatch between key and input sample names.
                            pearson_coeff = 'X'
                            spearman_coeff = 'X'

                        results[level].append((name, path.split('/')[-1],
                                                pearson_coeff, spearman_coeff))

    return results
