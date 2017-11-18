import os

from click.testing import CliRunner

from metacov.cli import pileup, scan

import pytest

import data

runner = CliRunner()


@pytest.mark.parametrize(
    'infile',
    [[data.fq1], [data.fq2], [data.fq1, data.fq2], [data.bam]]
)
def test_scan(tmpdir, infile):
    with tmpdir.as_cwd():
        result = runner.invoke(scan, infile + ["-o", "out1.csv", "-b", "out2.csv"])
        assert result.exit_code == 0, "exit code = "+str(result.exit_code)+" \n===\n" + result.output + "\n==="


def test_pileup(tmpdir):
    with tmpdir.as_cwd():
        result = runner.invoke(pileup, ["-b", data.bam,
                                        "-rb", data.regions_blast,
                                        "-o", "out.csv"])
        assert result.exit_code == 0, "\n===\n" + result.output + "\n==="


def test_pileup2(tmpdir):
    with tmpdir.as_cwd():
        result = runner.invoke(pileup, ["-b", data.bam,
                                        "-o", "out.csv"])
        assert result.exit_code == 0, "\n===\n" + result.output + "\n==="


def test_pileup_histogram(tmpdir):
    with tmpdir.as_cwd():
        result = runner.invoke(scan, [data.bam, "-o", "out1.csv"])
        assert result.exit_code == 0, "\n===\n" + result.output + "\n==="

        result = runner.invoke(pileup, ["-b", data.bam,
                                        "-rb", data.regions_blast,
                                        "-k", "out1.csv",
                                        "-o", "out2.csv"])
        assert result.exit_code == 0, "\n===\n" + result.output + "\n==="
