import os

from click.testing import CliRunner

from metacov.cli import pileup, scan

import pytest

runner = CliRunner()

path = os.path.join(os.path.dirname(__file__), 'data')
fq1 = os.path.join(path, 'ecoli_1K_1.fq.gz')
fq2 = os.path.join(path, 'ecoli_1K_2.fq.gz')
bam = os.path.join(path, 'bbmap.sorted.bam')
regions_blast = os.path.join(path, 'regions.blast7')

@pytest.mark.parametrize('infile', [[fq1], [fq2], [fq1, fq2], [bam]])
def test_scan(tmpdir, infile):
    with tmpdir.as_cwd():
        result = runner.invoke(scan, infile + ["-o", "out1.csv", "-b", "out2.csv"])
        assert result.exit_code == 0, "exit code = "+str(result.exit_code)+" \n===\n" + result.output + "\n==="


def test_pileup(tmpdir):
    with tmpdir.as_cwd():
        result = runner.invoke(pileup, ["-b", bam,
                                        "-rb", regions_blast,
                                        "-o", "out.csv"])
        assert result.exit_code == 0, "\n===\n" + result.output + "\n==="


def test_pileup2(tmpdir):
    with tmpdir.as_cwd():
        result = runner.invoke(pileup, ["-b", bam,
                                        "-o", "out.csv"])
        assert result.exit_code == 0, "\n===\n" + result.output + "\n==="


def test_pileup_histogram(tmpdir):
    with tmpdir.as_cwd():
        result = runner.invoke(scan, [bam, "-o", "out1.csv"])
        assert result.exit_code == 0, "\n===\n" + result.output + "\n==="

        result = runner.invoke(pileup, ["-b", bam,
                                        "-rb", regions_blast,
                                        "-k", "out1.csv",
                                        "-o", "out2.csv"])
        assert result.exit_code == 0, "\n===\n" + result.output + "\n==="
