import os
import sys


def require_file_exists_and_test_empty(infile):
    '''determines if an input file is present and exits with error message if
    untrue; returns whether the file is empty'''
    if os.path.exists(infile):
        if os.stat(infile).st_size == 0:
            return True
        return False
    sys.stderr.write('ERROR: {} file absent\n'.format(infile))
    sys.exit(1)
