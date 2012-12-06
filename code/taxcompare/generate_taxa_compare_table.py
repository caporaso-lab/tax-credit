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
from os.path import exists
from qiime.workflow import WorkflowError
from qiime.parse import parse_taxa_summary_table
from qiime.compare_taxa_summaries import compare_taxa_summaries

assignment_method_choices = ['rdp','blast','rtax','mothur','tax2tree']

def format_output(compare_tables, separator):
    """Formats the output from generate_taxa_compare_table into a list of lists(corresponding to files) ready for writelines()"""
    result = []
    for i, table in enumerate(compare_tables):
        result.append(list())
        datasets = sorted(table.keys())
        methods = set()
        for m in table.itervalues():#Find all methods used in table
            methods|= set(m.keys())
        methods = sorted(list(methods))
        result[i].append('P'+separator+'S\t'+'\t'.join(methods)+'\n')
        for dataset in datasets:
            line = dataset+'\t'
            for method in methods:
                try:
                    line += table[dataset][method][0] + separator + table[dataset][method][1]+'\t'
                except KeyError:
                    line += 'N/A'+'\t' #Don't have data for that set/method
            result[i].append(line + '\n')
    return result

def get_key_files(directory):
    """Given a directory containing keys, will identify key files and return a dict {name of study: file path}"""
    if(not exists(directory)):
        raise WorkflowError('The key directory does not exist.')
    key_fps = {}
    for key_file in walk(directory).next()[2]:
        study = key_file.split('_')[0].capitalize()
        key_fps[study] = directory+'/'+key_file
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

def generate_taxa_compare_table(root, key_directory, levels=[2,3,4,5,6]):
    """Walks a file tree rooted at root, finds otu tables and compares them against """\
    """the keys in key_directory. Will check every otu table at level in levels"""
    key_fps = get_key_files(key_directory)

    results = []

    
    for i in range(len(levels)):
        results.append(dict())

    for(path, dirs, files) in walk(root):
        for choice in assignment_method_choices:
            if choice in path:
                experiment = path.split('/')[-2].rstrip('-123').capitalize()
                for f in files:
                    if 'otu_table_mc2_w_taxa_L' in f:
                        name = path.split('/')[-2].capitalize()
                        level = int(f[-5])

                        with open(path+'/'+f,'U') as run_file:
                            #print run_file.readlines()
                            #run_file.seek(0)
                            test = run_file.readline()
                            if('Taxon\t' not in test):
                                raise WorkflowError('Invalid multiple_assign_taxonomy output file, check for corrupted file: '+path)
                            run_file.seek(0)
                            run = parse_taxa_summary_table(run_file)

                        with open(key_fps[experiment]) as key_file:
                            test = key_file.readline()
                            if('Taxon\t' not in test):
                                raise WorkflowError('Invalid key file in directory: '+path)
                            key_file.seek(0)
                            key = parse_taxa_summary_table(key_file)

                        try:
                            pearson_coeff, spearman_coeff = get_coefficients(run, key)
                        except ValueError:#compare_taxa_summaries couldn't find a match between the 2
                        #Likely due to mismatch between key and input sample names.
                            pearson_coeff = 'X'
                            spearman_coeff = 'X'
                        try:
                            results[level-2][name]
                        except KeyError:
                            results[level-2][name] = dict()
                        results[level-2][name][path.split('/')[-1]] = (pearson_coeff, spearman_coeff)
    return results
