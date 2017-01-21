#!/usr/bin/env python


# ----------------------------------------------------------------------------
# Copyright (c) 2016--, tax-credit development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from skbio import io
from re import search
from math import ceil
from collections import OrderedDict, Counter
from random import shuffle, choice
from os.path import isfile, exists, join, basename
from os import path, makedirs


# Sniff object: if file import to list, if list returns list, else error
def accept_list_or_file(infile, field=None, delim='\t'):
    '''if file return list, if list return list, else error'''
    if isinstance(infile, list):
        if field is not None:
            listout = [line.split(delim)[field].strip() for line in infile]
        else:
            listout = infile
    elif isfile(infile):
        listout = import_to_list(infile, field=field, delim=delim)
    return listout


# Import lines of text file into list object
def import_to_list(infile, field=None, delim='\t'):
    '''File -> list object'''

    with open(infile, "r") as inputfile:
        if field is None:
            lines = [line.strip() for line in inputfile]
        else:
            lines = [line.strip().split(delim)[field] for line in inputfile]
    return lines


# Import lines of text file into list object
def import_taxonomy_to_dict(infile):
    ''' taxonomy file -> dict'''
    with open(infile, "r") as inputfile:
        lines = {line.strip().split('\t')[0] : line.strip().split('\t')[1]\
                 for line in inputfile}
    return lines


# Write text file from dictionary
def export_list_to_file(input_list, outfile):
    '''list -> file'''

    with open(outfile, "w") as printout:
        printout.write('\n'.join(str(v) for v in input_list))


def output_list_or_file(listin, outfile):
    '''return list as list or file'''
    if outfile is not False:
        export_list_to_file(listin, outfile)
    else:
        return listin


def extract_rownames(infile):
    '''Extract seq ids (rownames) from file OR LIST'''
    # sniff filtermap and import to list if file
    line_list = accept_list_or_file(infile)

    # Extract rownames from list and pass to set
    id_list = {line.split('\t')[0].rstrip() for line in line_list}
    return id_list


# filter fasta based on taxonomy strings
def filter_sequences(infile, outfile, filtermap, keep=True):
    '''Filter a fasta file to contain (or exclude) sequnces matching search
    criteria listed in 'filtermap' (list of ids or tsv with ids as rownames)
    keep = Toggle whether to keep only sequnces listed in 'filtermap' or
           exclude these sequences from the output fasta.
    '''
    # Extract sequence IDs from filtermap
    id_list = extract_rownames(filtermap)

    # Loop infile, only keep (or exclude) sequences with matching IDs
    with open(outfile, 'w') as output_fasta:
        for seq in io.read(infile, format='fasta'):
            if seq.metadata['id'] in id_list and keep is True \
               or seq.metadata['id'] not in id_list and keep is False:
                seq.write(output_fasta, format='fasta')


# Gather seq IDs from fasta file headers
def extract_fasta_ids(infile):
    '''Gather sequence IDs from the header lines of a fasta file and output
    as a set. file -> set
    '''

    id_list = set()
    for sequence in io.read(infile, format='fasta'):
        id_list.add(sequence.metadata['id'])
    return id_list


# filter taxonomy strings on fasta (seq ID) or taxonomy (search terms)
def string_search(infile, pattern, discard=False, field=slice(None),
                  delim=';'):
    '''Search lines of file for pattern(s). Retain (default) or discard
    (discard = True) matching lines. Returns a new file containing matching /
    non-matching lines.

    <pattern> = one or more search patterns in |-separated string, enclosed in
                a single set of single quotes. E.g., to search 'a','b', or 'c',
                use: 'a|b|c'

    field = delim-delimited field to search. If field is not slice(None), will
        only search delim-delimited field. Otherwise will search whole line by
        default.

    FIELD MUST USE SLICE NOTATION or will FAIL. E.g., to search field 6, use:
        field = slice(6,7)

    infile = file or list object
    [file ->] list -> list
    '''

    # sniff filtermap and import to list if file
    search_list = accept_list_or_file(infile)

    # Keep or discard lines matching pattern
    keep_list = [line for line in search_list
                 if search(pattern, delim.join(line.split(delim)[field]))\
                 and discard is False \
                 or not search(pattern, delim.join(line.split(delim)[field]))\
                 and discard is True]

    return keep_list


# generate expected taxonomy files for novel taxa
def trim_taxonomy_strings(infile, level, delim=';'):
    '''Generate expected taxonomy strings for 'novel taxa'
    from input file tab-sep format (tsv) OR list of lines from tsv:
    seqID    kingdom;phylum;class;order;family;genus;species

    Outputs trimmed taxonomy strings up to specified level
    where 0 = kingdom, 1 = phylum, 2 = class, 3 = order, 4 = family,
    5 = genus, 6 = species.

    infile = file OR list
    filestem = filestem/basename without extension; will append '-L[level].tsv'
        if filestem = False, output defaults to nested list
    '''

    # sniff filtermap and import to list if file
    input_list = accept_list_or_file(infile)
    # Generate expected taxonomies by slicing to level
    taxa_list = [delim.join(line.split(delim)[0:level+1])
                 for line in input_list]
    return taxa_list


# sample unique taxonomy strings
def unique_lines(infile, mode='n', field=None, delimiter='\t',
                 printfield=False):
    '''Extract unique (taxonomy) strings from list or infile in format:
    seqID    kingdom;phylum;class;order;family;genus;species

    Match by 'delimiter'-separated field number 'field' (leftmost column = 0).
        field = None (default) will evaulate whole lines.

    mode:
    n = Output one copy of each unique line to outfile. (Default)
    u = Only output lines that are not repeated in the input.
    d = Only output lines that are repeated in the input.

    printfield: print only matching field (True) or whole line (False).
        If field = None, printfield will always print whole line.
    '''

    # sniff filtermap and import to list if file
    input_list = accept_list_or_file(infile)

    record = OrderedDict()
    record_counter = Counter()

    for line in input_list:
        printline = line
        # Toggle to match only 'delimiter'-sep column number 'field'.
        if field is not None:
            line = line.split(delimiter)[field]
            if printfield is True:
                printline = line

        # Add to record and counter
        record_counter.update({line: 1})
        if line not in record:
            record[line] = printline

    print_list = [record[key] for key, value in record.items() if mode == 'n' \
                  or (mode == 'd' and record_counter[key] > 1) \
                  or (mode == 'u' and record_counter[key] == 1)]

    return print_list


# remove non-branching taxa
def branching_taxa(infile, field=None, delim=';'):
    '''A modification of unique_lines(); this version detects branching
    taxonomic lineages. A major distinction is that this function prints ALL
    duplicates and no uniques, rather than a representative entry for
    duplicated entries (following unix uniq function).

    1) Read lines, store to memory as printline
        U33070	Root;Basidiomycota;Agaricomycetes;Agaricales
        U33090	Root;Basidiomycota;Agaricomycetes;Atheliales
        AF178417	Root;Ascomycota;Sordariomycetes;Hypocreales
    2) Cut off seq ID.
        Root;Basidiomycota;Agaricomycetes;Agaricales
        Root;Basidiomycota;Agaricomycetes;Agaricales
        Root;Ascomycota;Sordariomycetes;Hypocreales
    3) Cut off terminal taxon name to create search line.
        Root;Basidiomycota;Agaricomycetes
        Root;Basidiomycota;Agaricomycetes
        Root;Ascomycota;Sordariomycetes
    4) Print printline if search line is duplicated.
        U33070	Root;Basidiomycota;Agaricomycetes;Agaricales
        U33090	Root;Basidiomycota;Agaricomycetes;Atheliales

        NOT YET UNIT TESTED
    '''

    # sniff filtermap and import to list if file
    input_list = accept_list_or_file(infile)

    record = OrderedDict()
    record_counter = Counter()

    for printline in input_list:
        # Cut off seq ID
        searchline = printline.split('\t')[1]
        # Cut off terminal taxon name to create search line.
        if field is not None:
            searchline = delim.join(searchline.split(delim)[0:field])
        # push printline (always unique) to record
        record[printline] = searchline
        record_counter.update({searchline: 1})

    # Print if searchlines (values in record) are duplicated
    print_list = [key for key, value in record.items()
                  if record_counter[record[key]] > 1]

    return print_list


def stratify_taxonomy_subsets(infile, number_of_splits, basedir,
                              basename, level, delim=';'):
    '''[file] -> list -> files
    Split an input file or list into x parts, where x = number of splits.
    Do so in such a way that taxa at 'delim'-delimited field 'level' are
    evenly stratified across output files. Each part will not necessarily
    contain an even number of lines, though randomness should smooth this out.

    basedir = output directory that will contain output directories
    basename = name of directories that will be created in basedir in format
               basedir/basename-iter0/
    '''
    # sniff filtermap and import to list if file
    line_list = accept_list_or_file(infile)

    # input list is split into N chunks,
    chunk_choices = {i for i in range(number_of_splits)}
    chunks = [[] for i in range(number_of_splits)]

    # evenly distribute taxa at specified level among chunks
    taxon_dict = dict()
    for taxonomy in line_list:
        taxon = taxonomy.split(delim)[level]
        # assign taxon to random chunk
        chunk_choice = choice(list(chunk_choices))
        if taxon not in taxon_dict:
            taxon_dict[taxon] = set([chunk_choice])
        else:
            # if all chunk choices have not been called for taxon
            if len(taxon_dict[taxon]) < len(chunk_choices):
                # choose a chunk choice that has not been called
                chunk_choice = choice(list(chunk_choices - taxon_dict[taxon]))
            taxon_dict[taxon].add(chunk_choice)
        chunks[chunk_choice].append(taxonomy)

    # Generate output files, write N lines, where N=lines_per_file
    for i in range(0, len(chunks)):
        outdir = path.join(basedir, '{0}-iter{1}'.format(basename, i))
        if not exists(outdir):
            makedirs(outdir)
        export_list_to_file(chunks[i], path.join(outdir, 'query_taxa.tsv'))


def file_splitter(infile, number_of_splits, basedir, basename, random=True):
    '''Split an input file OR LIST into x parts, where x = number_of_splits.
    If random=True, lines will be randomly split to output files.
    If random=False, lines will be split into x contiguous chunks.

    basedir = output directory that will contain output directories
    basename = name of directories that will be created in basedir in format
               basedir/basename-iter0/

    NOT UNIT TESTED
    '''

    # sniff filtermap and import to list if file
    line_list = accept_list_or_file(infile)

    # Define number of lines to write to each split
    lines_per_file = ceil(len(line_list) / number_of_splits)

    if random is True:
        shuffle(line_list)

    # Break list into chunks
    file_chunks = [line_list[i:i + lines_per_file]
                   for i in range(0, len(line_list), lines_per_file)]

    # Generate output files, write N lines, where N=lines_per_file
    for i in range(0, len(file_chunks)):
        outdir = path.join(basedir, '{0}-iter{1}'.format(basename, i))
        if not exists(outdir):
            makedirs(outdir)
        export_list_to_file(file_chunks[i], path.join(outdir, 'query_taxa.tsv'))


def extract_taxa_names(infile, level=slice(6, 7), l_delim=';', field=None,
                       f_delim='\t'):
    '''Extract taxon names at a given level from taxonomy file OR LIST
    field = taxonomic level of interest. 0 = kingdom, 1 = phylum, 2 = class,
    3 = order, 4 = family, 5 = genus, 6 = species.
    Must use slice notation. For species, use level=slice(6, 7)
    NOT UNIT TESTED
    '''
    # sniff filtermap and import to list if file
    line_list = accept_list_or_file(infile, field=field, delim=f_delim)

    # Truncate taxonomies and pass to set
    name_list = [l_delim.join(line.split(l_delim)[level]).strip()\
                 for line in line_list]
    return name_list
