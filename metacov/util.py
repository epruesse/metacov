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

    colnums = []
    for names in (('sacc', 'sequence_id'),
                  ('sstart', 'start'),
                  ('send', 'end', 'stop')):
        col = None
        for name in names:
            if name in columns:
                col = columns.index(name)
                break
        if col is None:
            raise ValueError("Region file must have a column with a name in {}"
                             "".format(names))
        colnums.append(col)

    for row in csv_reader:
        yield Region("", row[colnums[0]], row[colnums[1]], row[colnums[2]])


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
