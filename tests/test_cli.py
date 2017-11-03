import os

from click.testing import CliRunner

from metacov.cli import pileup, scan

import pytest

runner = CliRunner()

path = os.path.join(os.path.dirname(__file__), 'data')
fq1 = os.path.join(path, 'ecoli_1K_1.fq.gz')
fq2 = os.path.join(path, 'ecoli_1K_2.fq.gz')


def test_scan(tmpdir):
    with tmpdir.as_cwd():
        result = runner.invoke(scan, [fq1, "out"])
        assert result.exit_code == 0, "\n===\n" + result.output + "\n==="
