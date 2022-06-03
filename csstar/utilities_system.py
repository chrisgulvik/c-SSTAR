#!/usr/bin/env python

import logging
import os
import subprocess
import sys


def check_python_version():
    if sys.version_info.major == 2 and sys.version_info.minor < 7:
        sys.stderr.write('Error: Python 2.7 or later required\n')
        sys.exit(1)
    elif sys.version_info.major == 3 and sys.version_info.minor < 3:
        sys.stderr.write('Error: Python 3.3 or later required\n')
        sys.exit(1)


def require_dependency(dep):
    for path in os.environ.get('PATH', '').split(':'):
        dep_path = os.path.join(path, dep)
        if os.path.exists(dep_path) and not os.path.isdir(dep_path):
            return True
    sys.stderr.write('ERROR: {} unavailable; not in $PATH\n'.format(dep))
    sys.exit(1)


def sys_call(syscmd):
    '''runs a system command string in a subshell and exits with an error
    message if the return code of the command was anything other than 0'''
    with open(os.devnull) as dump:
        returncode = subprocess.call(syscmd, stdout=dump, stderr=dump,
                                     shell=True)
        if returncode != 0:
            logging.error('failed sys_call ' + syscmd)
            sys.stderr.write('ERROR: failed sys_call {}\n'.format(syscmd))
            sys.exit(1)
