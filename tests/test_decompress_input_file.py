#!/usr/bin/env python

import os
from shutil import rmtree
from tempfile import mkdtemp

from csstar.file_convert import decompress_file


def test_decompress_input_file(capsys):
    tmp = mkdtemp()
    infile = os.path.join(
        os.path.dirname(__file__), 'data', 'GCF_013463375.1.faa.gz')
    faa = decompress_file(infile, tmp)
    assert os.path.isfile(faa)
    assert os.path.getsize(faa) > 0
    with open(faa) as ifh:
        first_line = ifh.readline()
        for line in ifh:
            pass
        last_line = line.rstrip()
    assert first_line.startswith('>WP_020915717.1')
    assert last_line.endswith('INNKINCKTIAGGIGIDRLFFYIKKLKTIKKVYD')
    rmtree(tmp)
