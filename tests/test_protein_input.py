#!/usr/bin/env python

import os
import pytest
from shutil import rmtree
from tempfile import mkdtemp

from csstar.file_content import evaluate_nucleotide_composition
from csstar.file_convert import decompress_file


def input_file_composition(infile, min_ACGT, capsys, expect=0):
    with pytest.raises(SystemExit) as e:
        evaluate_nucleotide_composition(infile, min_ACGT)
    std_out, std_err = capsys.readouterr()
    assert std_err.startswith('ERROR:')
    assert e.type == SystemExit
    assert e.value.code == expect


def test_protein_input(capsys):
    tmp = mkdtemp()
    infile = os.path.join(
        os.path.dirname(__file__), 'data', 'GCF_013463375.1.faa.gz')
    faa = decompress_file(infile, tmp)
    input_file_composition(faa, 97, capsys, expect=1)
    rmtree(tmp)
