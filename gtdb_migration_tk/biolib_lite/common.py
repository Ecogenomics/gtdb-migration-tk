###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

"""
Those functions are copied from `biolib <https://github.com/dparks1134/biolib>`_  on the 09/10/2020

"""

import os
import errno
import sys
import logging
import ntpath
import mmap
import re
import gzip
from datetime import datetime

from numpy import (abs as np_abs,
                   array as np_array)





def canonical_gid(gid):
    """Get canonical form of NCBI genome accession.

    Example:
        G005435135 -> G005435135
        GCF_005435135.1 -> G005435135
        GCF_005435135.1_ASM543513v1_genomic -> G005435135
        RS_GCF_005435135.1 -> G005435135
        GB_GCA_005435135.1 -> G005435135
    """

    if gid.startswith('U'):
        return gid

    gid = gid.replace('RS_', '').replace('GB_', '')
    gid = gid.replace('GCA_', 'G').replace('GCF_', 'G')
    if '.' in gid:
        gid = gid[0:gid.find('.')]

    return gid


def get_num_lines(file_path):
    """ Calculate the number of lines in a file."""
    fp = open(file_path, "r+")
    buf = mmap.mmap(fp.fileno(), 0)
    lines = 0
    while buf.readline():
        lines += 1
    return lines


def is_float(s):
    """Check if a string can be converted to a float.
    Parameters
    ----------
    s : str
        String to evaluate.
    Returns
    -------
    boolean
        True if string can be converted, else False.
    """

    try:
        float(s)
    except ValueError:
        return False

    return True


def check_file_exists(input_file):
    """Check if file exists."""
    if not os.path.exists(input_file) or not os.path.isfile(input_file):
        logger = logging.getLogger('timestamp')
        logger.error('Input file does not exists: ' + input_file + '\n')
        sys.exit()


def check_dir_exists(input_dir):
    """Check if directory exists."""
    if not os.path.exists(input_dir) or not os.path.isdir(input_dir):
        logger = logging.getLogger('timestamp')
        logger.error('Input directory does not exists: ' + input_dir + '\n')
        sys.exit()


def make_sure_path_exists(path):
    """Create directory if it does not exist."""

    if not path:
        # lack of a path qualifier is acceptable as this
        # simply specifies the current directory
        return

    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            logger = logging.getLogger('timestamp')
            logger.error('Specified path could not be created: ' + path + '\n')
            sys.exit()


def remove_extension(filename, extension=None):
    """Remove extension from filename.
    A specific extension can be specified, otherwise
    the extension is taken as all characters after the
    last period.
    """
    f = ntpath.basename(filename)

    if extension and f.endswith(extension):
        f = f[0:f.rfind(extension)]
    else:
        f = os.path.splitext(f)[0]

    if f[-1] == '.':
        f = f[0:-1]

    return f

#### PC ADDED functionalities to common #####

def select_delimiter(metafile):
    # Parse TSV or CSV file
    for line in open(metafile):
        if len(line.split('\t')) >= len(line.split(',')):
            return '\t'
        else:
            return ','

def clean_html(raw_html):
    """
    Remove all HTML tags from HTML line
    """
    cleanr = re.compile('<.*?>')
    cleantext = re.sub(cleanr, '', raw_html)
    return cleantext.strip()

def clean_parenthesis(raw_line):
    """
    Remove parenthesis content from HTML line
    """
    # result = re.sub("[\(\[].*?[\)\]]", "", raw_line)
    result = re.sub("\(see[^\)]+\)", "", raw_line)
    return result
