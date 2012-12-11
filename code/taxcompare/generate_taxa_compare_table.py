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
from os.path import exists, join
from qiime.workflow import WorkflowError
from qiime.parse import parse_taxa_summary_table
from qiime.compare_taxa_summaries import compare_taxa_summaries

assignment_method_choices = ['rdp','blast','rtax','mothur','tax2tree']

def format_output(compare_tables, separator):
    """Formats the output from generate_taxa_compare_table into a list of lists(corresponding to files) ready for writelines()"""
    result = []
    for i, table in enumerate(compare_tables):
        result.append(list())
        if not table:
            continue
        datasets = sorted(table.keys())
        methods = set()
        for m in table.itervalues():
            #Find all methods used in table
            methods|= set(m.keys())
        methods = sorted(list(methods))
        result[i].append('P'+separator+'S\t'+'\t'.join(methods)+'\n')
        for dataset in datasets:
            line = dataset+'\t'
            for method in methods:
                try:
                    line += table[dataset][method][0] + separator + table[dataset][method][1]+'\t'
                except KeyError:
                    #Don't have data for that set/method
                    line += 'N/A'+'\t'
            result[i].append(line + '\n')
    return result

def get_key_files(directory):
    """Given a directory containing keys, will identify key files and return a dict {name of study: file path}"""
    if(not exists(directory)):
        raise WorkflowError('The key directory does not exist.')
    key_fps = {}
    for key_file in walk(directory).next()[2]:
        if key_file.endswith('~'):
            continue
        study = key_file.split('_')[0].capitalize()
        key_fps[study] = join(directory, key_file)
    if(not key_fps):
        raise WorkflowError('There are no key files in the given directory.')
    return key_fps

def get_coefficients(run, key):
    """Given a parsed taxa summary table, will find and return correlation coefficients"""
    pearson_compare = compare_taxa_summaries(run, key, 'paired', 'pearson')
    spearman_compare = compare_taxa_summaries(run, key, 'paired', 'spearman')

    pearson_coeff = pearson_compare[2].split('\n')[-2].split()[0]
    spearman_coeff = spearman_compare[2].split('\n')[-2].split()[0]

    return pearson_coeff, spearman_coeff

def generate_taxa_compare_table(root, key_directory, levels=None):
    """Finds otu tables in root and compares them against the keys in key_directory.

    Walks a file tree starting at root and finds the otu tables output by
    multiple_assign_taxonomy.py. Then compares the found otu tables to their corresponding
    key in key_directory. Returns a list containing a dict for every
    level of output compared. It should be noted that since levels start at 2, by default 
    a specific level will be at results[level-2]. Each level dict is of the format 
    {name of study: {method_and_params: (pearson, spearman)}}

    Parameters:
    root: path to root of multiple_assign_taxonomy.py output.
    key_directory: path to directory containing known/expected compositions. Each study
        should be in its own otu table.
    levels: INCOMPLETE. Use other than default will cause unexpected results. The
        multiple_assign_taxonomy.py output levels to be analyzed."""
    key_fps = get_key_files(key_directory)

    results = []

    if not levels:
        levels = [2,3,4,5,6]

    if len(levels) > 5:
        raise WorkflowError('Too many levels.')
    for l in levels:
        if l < 2 or l > 6:
            raise WorkflowError('Level out of range: ' + str(l))

    for i in range(5):
        results.append(dict())

    for(path, dirs, files) in walk(root):
        for choice in assignment_method_choices:
            #Checks if this dir's name includes a known assignment method (and therefor contains that output)
            if choice in path:
                study = path.split('/')[-2].rstrip('-123').capitalize()
                for f in files:
                    if 'otu_table_mc2_w_taxa_L' in f and not f.endswith('~'):
                        name = path.split('/')[-2].capitalize()
                        level = int(f[-5])
                        if level not in levels:
                            #If that level wasn't requested, skip it.
                            continue

                        with open(join(path,f),'U') as run_file:
                            #Open and parse run file
                            test = run_file.readline()
                            if('Taxon\t' not in test):
                                raise WorkflowError('Invalid multiple_assign_taxonomy output file, check for corrupted file: '+path)
                            run_file.seek(0)
                            run = parse_taxa_summary_table(run_file)

                        with open(key_fps[study],'U') as key_file:
                            #Open and parse key file
                            test = key_file.readline()
                            if('Taxon\t' not in test):
                                raise WorkflowError('Invalid key file in directory: '+path)
                            key_file.seek(0)
                            key = parse_taxa_summary_table(key_file)

                        try:
                            pearson_coeff, spearman_coeff = get_coefficients(run, key)
                        except ValueError:
                            #compare_taxa_summaries couldn't find a match between the 2
                            #Likely due to mismatch between key and input sample names.
                            pearson_coeff = 'X'
                            spearman_coeff = 'X'
                        try:
                            results[level-2][name]
                        except KeyError:
                            results[level-2][name] = dict()
                        results[level-2][name][path.split('/')[-1]] = (pearson_coeff, spearman_coeff)
    return results
