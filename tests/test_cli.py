import os

from click.testing import CliRunner

from metacov.cli import pileup, scan

import pytest

import data

runner = CliRunner()

def check_click_result(result):
    from traceback import format_exception
    if result.exception or result.exit_code != 0:
        print("Output: ", result.output)
    assert not result.exception, format_exception(*result.exc_info)
    assert result.exit_code == 0



@pytest.mark.parametrize(
    'infile',
    [[data.fq1], [data.fq2], [data.fq1, data.fq2], [data.bam]]
)
def test_scan(tmpdir, infile):
    with tmpdir.as_cwd():
        result = runner.invoke(scan, infile + ["-o", "out1.csv", "-b", "out2.csv"])
        check_click_result(result)


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
        check_click_result(result)
