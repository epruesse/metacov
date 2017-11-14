import click
import logging
import pysam
import csv
import os
from io import open

from metacov import pileup as _pileup
from metacov import scan as _scan
from metacov import blast
from metacov.pyfq import FastQFile, FastQFilePair


logging.basicConfig(
    level=logging.INFO,
    format="[%(relativeCreated)6.1f %(funcName)s]  %(message)s",
    datefmt="%I:%M:%S"
)

log = logging.getLogger(__name__)


@click.group()
def main():
    """
    Metacov -- Calculate abundance estimates over metagenome regions
    """


@main.command()
@click.argument('bamfile', type=click.File('rb'))
@click.argument('regionfile_name')  #, type=click.File('r'))
@click.argument('coveragefile', type=click.File('w'))
@click.option(
    '--kmer-histogram', '-k', type=click.File('r')
)
def pileup(bamfile, regionfile_name, coveragefile, kmer_histogram):
    """
    Compute coverage depths

    \b
    Arguments:
      BAMFILE      Input BAM file. Must be sorted and indexed
      REGIONFILE  Input Region file. Can be one of: BLAST 7
    """
    regionfile = open(regionfile_name, "r")
    log.info("Gathering coverages for regions in \"{}\""
             "".format(regionfile.name))
    blastreader = blast.reader(regionfile)

    log.info("Processing BAM file \"{}\"".format(bamfile.name))
    bam = pysam.AlignmentFile(bamfile.name)

    try:
        log.info("Number of reads:\n"
                 "  total:    {total}\n"
                 "  mapped:   {mapped} ({mappct}%)\n"
                 "  unmapped: {unmapped}\n"
                 "".format(
                     mapped=bam.mapped, unmapped=bam.unmapped,
                     total=bam.mapped+bam.unmapped,
                     mappct=bam.mapped / (bam.mapped+bam.unmapped) * 100
                 ))
    except AttributeError:
        log.error("BAM file not indexed!?")

    name2ref = {word.split()[0]: word
                for word in bam.references}

    k_cor = None
    if kmer_histogram is not None:
        k_cor = _pileup.load_kmerhist(kmer_histogram)

    writer = None

    regionfile_size = os.path.getsize(regionfile.name)
    with click.progressbar(
            # length=os.path.getsize(regionfile.name),
            length=0,
            label="Calculating coverages") as bar:
        for hit in blastreader:
            # translate hit name, in case the bam file contains
            # multi-word references
            ref = name2ref[hit.sacc]

            # sort the start/end (blast uses end<start for reverse strand)
            start, end = sorted((hit.sstart, hit.send))
            length = end-start

            result = _pileup.classic(bam, ref, start, end)

            if k_cor is not None:
                result.update(_pileup.experimental(bam, k_cor,
                                                   ref, start, end))

            if writer is None:
                fieldnames = ['sacc', 'start', 'end'] + sorted(list(result.keys()))
                writer = csv.DictWriter(coveragefile, fieldnames=fieldnames)
                writer.writeheader()

            result.update({
                'sacc': hit.sacc,
                'start': hit.sstart,
                'end': hit.send
            })
            writer.writerow(result)
            # bar.update(regionfile.tell())


@main.command()
@click.argument(
    "readfile", nargs=-1, required=True
)  # , type=click.File("r"))
@click.option(
    "--readfile-type", "-t", type=click.Choice(['bam', 'fq']),
    help=""
)
@click.option(
    "--out-basehist", "-b", type=click.File("w", lazy=False),
    help=""
)
@click.option(
    "--out-kmerhist", "-o", type=click.File("w", lazy=False)
)
@click.option(
    "-k", type=int, default=7,
    help="Length of k-mer"
)
@click.option(
    "--step", "-s", type=int, default=7,
    help=""
)
@click.option(
    "--offset", "-O", type=int, default=0,
    help=""
)
@click.option(
    "--number", "-n", type=int, default=8,
    help=""
)
def scan(readfile, readfile_type, out_basehist, out_kmerhist,
         k, number, step, offset):
    # Try to guess readfile type if none specified
    if not readfile_type:
        for ext, ft in {'.bam': 'bam',
                        '.sam': 'bam',
                        '.fq': 'fq',
                        '.fq.gz': 'fq',
                        '.fastq': 'fq',
                        '.fastq.gz': 'fq'}.items():
            if readfile[0].endswith(ext):
                readfile_type = ft
                break

    # Fail if readfile not recognized above
    if not readfile_type:
        raise click.UsageError(
            "Couldn't guess input format. Please supply -t"
        )

    # Check only one bam file
    if readfile_type != 'fq' and len(readfile) > 1:
        raise click.UsageError(
            "Multiple input files only supported for fastq"
        )

    # Check at most two read files
    if len(readfile) > 2:
        raise click.UsageError(
            "At most two fastq files allowed (fwd and rev)"
        )

    update_every = 100000
    max_reads = update_every * 10
    max_reads = 0

    if readfile_type == 'bam':
        infile = pysam.AlignmentFile(readfile[0])
        log.info("mapped = {}, unmapped = {}, total = {}"
                 "".format(
                     infile.mapped, infile.unmapped,
                     infile.mapped + infile.unmapped))
        length = infile.mapped + infile.unmapped

        def update_func(x):
            return x.update(update_every)

    elif readfile_type == 'fq':
        if len(readfile) > 1:
            infile = FastQFilePair(readfile[0], readfile[1])
        else:
            infile = FastQFile(readfile[0])
        length = infile.size
        lastpos = 0

        def update_bar(x):
            nonlocal lastpos
            nonlocal infile
            x.update(infile.pos - lastpos)
            lastpos = infile.pos
        update_func = update_bar
    else:
        raise Exception("Internal error: Unknown readfile type")

    counters = []
    if out_basehist:
        counters.append(_scan.BaseHist())
    if out_kmerhist:
        counters.append(_scan.KmerHist(k, number, step, offset))

    counters = _scan.ByFlag(counters, [_scan.FLAG_READDIR, _scan.FLAG_MAPPED])

    with click.progressbar(length=length,
                           label="Scanning reads",
                           show_pos=True) as bar:
        with infile:
            nreads = _scan.scan_reads(
                infile, counters,
                update_every, lambda: update_func(bar),
                max_reads)
        update_func(bar)

    log.info("Processed {} reads".format(nreads))

    n = 0
    if out_basehist:
        out = csv.writer(out_basehist)
        out.writerows(counters.get_rows(n))
        n = n + 1

    if out_kmerhist:
        out = csv.writer(out_kmerhist)
        out.writerows(counters.get_rows(n))
        n = n + 1
