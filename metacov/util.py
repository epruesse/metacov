import csv
import datetime
from collections import namedtuple

import click

from metacov import blast


class Timeit:
    def __init__(self, message="timed: {elapsed}s"):
        self.message = message

    def __enter__(self):
        self.begin = datetime.datetime.now()
        return self

    def __exit__(self, type, value, traceback):
        self.stamp()

    def stamp(self):
        print(self.message.format(
            elapsed=datetime.datetime.now() - self.begin
        ))


Region = namedtuple("Region", ["qacc", "sacc", "sstart", "send"])


def get_regions_from_blast7(regionfile):
    """Yields Region tuples from BLAST file"""
    for hit in blast.reader(regionfile):
        yield hit


def get_regions_from_csv(regionfile):
    """Yields Region tuples from CSV file

    The column containing the contig name must be 'sacc' or 'sequence_id'.
    The column containing the start offset must be 'sstart' or 'start'.
    The column containing the end offset must be 'end' or 'send'.
    """
    csv_reader = csv.reader(regionfile)
    columns = next(csv_reader)

    # rename columns if needed, aliasing
    # sequence_id -> sacc
    # start -> sstart
    # end -> send
    if 'sacc' not in columns:
        if 'sequence_id' in columns:
            columns[columns.index('sequence_id')] = 'sacc'
        else:
            raise ValueError("can't find sequence_id/sacc in csv")

    if 'sstart' not in columns:
        if 'start' in columns:
            columns[columns.index('start')] = 'sstart'
        else:
            raise ValueError("can't find sstart/start in csv")

    if 'send' not in columns:
        if 'end' in columns:
            columns[columns.index('end')] = 'send'
        else:
            raise ValueError("can't find send/end in csv")

    for line in csv_reader:
        yield Region(line)


def get_regions_from_bam(bamfile):
    """Yields Region tuples directly from BAM file
    """
    for n, (rlen, rname) in enumerate(zip(bamfile.lengths,
                                          bamfile.references)):
        yield Region(n, rname, 0, rlen)


def make_region_iterator(regionfile_blast7, regionfile_csv, bam):
    # check params
    if regionfile_blast7 and regionfile_csv:
        raise click.BadParameter(
            "Only one of regionfile-blast7 and regionfile-csv may be specified"
        )

    if regionfile_blast7:
        return get_regions_from_blast7(regionfile_blast7)
    if regionfile_csv:
        return get_regions_from_csv(regionfile_csv)
    return get_regions_from_bam(bam)
