import os

from metacov.pyfq import FastQFile, FastQFilePair, FastQWriter

from pysam import AlignmentFile

import pytest

import data


@pytest.mark.parametrize(
    "fqfile,exp_size,exp_lengths,exp_base_freq",
    [
        (FastQFile(data.fq1),
         417,
         [0, 0, 0, 47, 88, 102, 134, 134, 231, 320, 1000],
         [44443, 45470, 44663, 43795, 2056]
         ),
        (FastQFile(data.fq2),
         414,
         [0, 0, 0, 69, 91, 104, 146, 144, 237, 318, 947],
         [44313, 44983, 43960, 42643, 2056]
         ),
        (FastQFilePair(data.fq1, data.fq2),
         831,
         [0, 0, 0, 116, 179, 206, 280, 278, 468, 638, 1947],
         [88756, 90453, 88623, 86438, 4112]
         ),
    ])
def test_FastQFile(fqfile, exp_size, exp_lengths, exp_base_freq):
    with fqfile as infile:
        assert infile.size == exp_size
        assert infile.pos == 0
        lengths = [0]*11
        base_freq = [0] * 5
        for read in infile:
            lengths[int(read.rlen/10)] += 1
            assert read.pos <= infile.size
            for base in read.seq:
                assert base < 5
                assert base >= 0
                base_freq[base] += 1
        assert exp_base_freq == base_freq
        assert exp_lengths == lengths


def test_FastQWriter(tmpdir):
    with tmpdir.as_cwd():
        n1 = 0
        with FastQFile(data.fq1) as infile, \
             FastQWriter("out.fq.gz") as outfile:
            for read in infile:
                outfile.write(read)
                n1 = n1 +1
        n2 = 0
        with FastQFile("out.fq.gz") as infile1, \
             FastQFile(data.fq1) as infile2:
            for read1, read2 in zip(infile1, infile2):
                assert read1.rlen == read2.rlen
                assert read1.seq == read2.seq
                n2 = n2 + 1
        assert n1 == n2, "n1 != n2 ({} != {})".format(n1, n2)
        assert n2 == 2056


def notest_FastQWriter2(tmpdir):
    fqfile = AlignmentFile(data.bam)
    fqwriter = FastQWriter("out.fq.gz")
    with tmpdir.as_cwd():
        n1 = 0
        with AlignmentFile(data.bam) as infile, \
             fqwriter as outfile:
            for read in infile:
                outfile.write(read)
                n1 = n1 +1
        n2 = 0
        with AlignmentFile(data.bam) as infile1, \
             FastQFile("out.fq.gz") as infile2:
            for read1, read2 in zip(infile1, infile2):
                assert read1.rlen == read2.rlen
                assert read1.seq == read2.seq
                print(read1.seq)
                n2 = n2 + 1
        assert n1 == n2, "n1 != n2 ({} != {})".format(n1, n2)
