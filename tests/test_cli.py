import os

from click.testing import CliRunner

from metacov.cli import pileup, scan

import pytest

runner = CliRunner()

path = os.path.join(os.path.dirname(__file__), 'data')
fq1 = os.path.join(path, 'ecoli_1K_1.fq.gz')
fq2 = os.path.join(path, 'ecoli_1K_2.fq.gz')
bam = os.path.join(path, 'bbmap.sorted.bam')

@pytest.mark.parametrize('infile', [[fq1], [fq2], [fq1, fq2], [bam]])
def test_scan(tmpdir, infile):
    with tmpdir.as_cwd():
        result = runner.invoke(scan, infile + ["out"])
        assert result.exit_code == 0, "\n===\n" + result.output + "\n==="
