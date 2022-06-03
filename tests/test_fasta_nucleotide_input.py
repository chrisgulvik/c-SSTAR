#!/usr/bin/env python

import os
import pytest
from shutil import rmtree
from tempfile import mkdtemp

from csstar.file_content import evaluate_nucleotide_composition
from csstar.file_convert import decompress_file


def test_fasta_nucleotide_input(capsys):
    tmp = mkdtemp()
    infile = os.path.join(
        os.path.dirname(__file__), 'data', 'SRR3112344.fna.gz')
    fna = decompress_file(infile, tmp)
    evaluate_nucleotide_composition(fna, 97)
    std_out, std_err = capsys.readouterr()
    assert len(std_err) == 0
    rmtree(tmp)
