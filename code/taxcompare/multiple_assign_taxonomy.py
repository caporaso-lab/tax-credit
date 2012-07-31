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

"""Contains functions used in the multiple_assign_taxonomy.py script."""

from qiime.workflow import (print_commands, print_to_stdout, no_status_updates,
                            call_commands_serially)

def assign_taxonomy_multiple_times(input_dir, assignment_methods,
        reference_seqs_fp, rep_set_filename, clean_otu_table_filename,
        rdp_id_to_taxonomy_fp=None, blast_id_to_taxonomy_fp=None,
        confidences=None, e_values=None,
        command_handler=call_commands_serially,
        status_update_callback=print_to_stdout, force=False):
    pass
