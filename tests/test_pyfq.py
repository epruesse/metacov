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
