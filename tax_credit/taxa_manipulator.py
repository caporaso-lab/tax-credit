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
from collections import OrderedDict, Counter
from random import choice
from os.path import isfile, exists
from os import path, makedirs
from biom import load_table
from biom.cli.util import write_biom_table


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
    else:
        raise RuntimeError('unknown file type in ' + infile)
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
        lines = {line.strip().split('\t')[0]: line.strip().split('\t')[1]
                 for line in inputfile}
    return lines


# Write text file from dictionary
def export_list_to_file(input_list, outfile):
    '''list -> file'''
    with open(outfile, "w") as printout:
        printout.write('\n'.join(str(v) for v in input_list))


def extract_rownames(infile):
    '''Extract seq ids (rownames) from file OR LIST'''
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
            if (seq.metadata['id'] in id_list and keep is True or
                    seq.metadata['id'] not in id_list and keep is False):
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
                  delim=';', f_field=None, f_delim='\t'):
    '''Search lines of file for pattern(s). Retain (default) or discard
    (discard = True) matching lines. Returns a new file containing matching /
    non-matching lines.

    <pattern> = one or more search patterns in |-separated string, enclosed in
                a single set of single quotes. E.g., to search 'a','b', or 'c',
                use: 'a|b|c'

    field = delim-delimited field to search. If field is not slice(None), will
        only search delim-delimited field. Otherwise will search whole line by
        default.

    f_field: f_delim-delimited field to search. If field is not None (default),
        whole lines will be extracted from infile, prior to delimited/searching
        with field and delim.

    FIELD MUST USE SLICE NOTATION or will FAIL. E.g., to search field 6, use:
        field = slice(6,7)

    infile = file or list object
    [file ->] list -> list
    '''
    search_list = accept_list_or_file(infile, field=f_field, delim=f_delim)

    # Keep or discard lines matching pattern
    keep_list = [line for line in search_list
                 if search(pattern, delim.join(line.split(delim)[field])) and
                 discard is False or
                 not search(pattern, delim.join(line.split(delim)[field])) and
                 discard is True]

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
    filestem = filestem/base without extension; will append '-L[level].tsv'
        if filestem = False, output defaults to nested list
    '''
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

    print_list = [record[key] for key, value in record.items()
                  if mode == 'n' or
                  (mode == 'd' and record_counter[key] > 1) or
                  (mode == 'u' and record_counter[key] == 1)]

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
    '''
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
                              base, level, delim=';'):
    '''[file] -> list -> files
    Split an input file or list into x parts, where x = number of splits.
    Do so in such a way that taxa at 'delim'-delimited field 'level' are
    evenly stratified across output files. Each part will not necessarily
    contain an even number of lines, though randomness should smooth this out.

    basedir = output directory that will contain output directories
    base = name of directories that will be created in basedir in format
               basedir/base-iter0/
    '''
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
        outdir = path.join(basedir, '{0}-iter{1}'.format(base, i))
        if not exists(outdir):
            makedirs(outdir)
        export_list_to_file(chunks[i], path.join(outdir, 'query_taxa.tsv'))


def extract_taxa_names(infile, level=slice(6, 7), l_delim=';', field=None,
                       f_delim='\t', stripchars=None):
    '''Extract taxon names at a given level from taxonomy file OR LIST
    field = taxonomic level of interest. 0 = kingdom, 1 = phylum, 2 = class,
    3 = order, 4 = family, 5 = genus, 6 = species.
    Must use slice notation. For species, use level=slice(6, 7)

    stripchars: str
        Default, None, will strip leading and trailing whitespace. Set to ""
        to turn off character stripping.
    '''
    line_list = accept_list_or_file(infile, field=field, delim=f_delim)

    # Truncate taxonomies and pass to set
    name_list = [l_delim.join(line.split(l_delim)[level]).strip(stripchars)
                 for line in line_list]
    return name_list


def compile_reference_taxa(ref_taxa_fp, delim='\t'):
    '''Extract taxa names as dict {taxon: id}. Accepts list or path as input'''
    ref_taxa = {}
    ref = accept_list_or_file(ref_taxa_fp)
    for l in ref:
        i, t = l.strip().split(delim)
        if t not in ref_taxa.keys():
            ref_taxa[t] = [i]
        else:
            ref_taxa[t].append(i)
    return ref_taxa


def convert_tsv_to_biom(infile, outfile, transpose=True,
                        obs_ids_to_metadata=False):
    table = load_table(infile)
    if transpose is True:
        table = table.transpose()
    if obs_ids_to_metadata is True:
        metadata = {sid: {'taxonomy': sid.split(';')}
                    for sid in table.ids(axis='observation')}
        table.add_metadata(metadata, 'observation')
    write_biom_table(table, 'hdf5', outfile)
