import csv
import logging
from collections import OrderedDict

import click

from metacov import pileup as _pileup
from metacov import scan as _scan
from metacov import util
from metacov.pyfq import FastQFile, FastQFilePair, FastQWriter

import numpy as np

import pysam

import sys

logging.basicConfig(
    level=logging.INFO,
    format="[%(relativeCreated)6.1f %(funcName)s]  %(message)s",
    datefmt="%I:%M:%S"
)

log = logging.getLogger(__name__)


@click.group()
def main():
    """
    MetaCov estimates abundance values from the stacking depth of
    reads mapped to a reference.
    """


@main.command()
@click.option('--bamfile', '-b', type=click.File('rb'), required=True,
              help="Input BAM file. Must be sorted and indexed.")
@click.option('--reference-fasta', '-f', type=click.File('rb'))
@click.option('--regionfile-blast7', '-rb', type=click.File('r'),
              help="Input Region file in BLAST7 format")
@click.option('--regionfile-csv', '-rc', type=click.File('r'),
              help="Input Region file in CSV format")
@click.option('--kmer-histogram', '-k', type=click.File('r'),
              help="Kmer Histogram produced with metacov scan")
@click.option('--kmer-length', '-K', type=int, default=7,
              help="Length of k-mer")
@click.option('--outfile', '-o', type=click.File('w'), default="-",
              help="Output CSV (default STDOUT)")
def pileup(bamfile, reference_fasta, regionfile_blast7, regionfile_csv,
           kmer_histogram, kmer_length, outfile):
    """
    Compute fold coverage values
    """

    # open bamfile
    bam = pysam.AlignmentFile(bamfile.name)

    # open reference
    fasta = pysam.FastaFile(reference_fasta.name) if reference_fasta else None
    # FIXME: If we get an OSError-fail-to-open here, it's because the file
    #        isn't formatted to spec (varying line length)

    # get region iterator
    regions = util.make_region_iterator(regionfile_blast7, regionfile_csv, bam)

    # dump stats
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

    name2ref = {word.split()[0]: word for word in bam.references}
    k_cor = _pileup.load_kmerhist(kmer_histogram) if kmer_histogram else None
    writer = None

    with click.progressbar(file=sys.stderr, length=0, label="Calculating coverages"):  # as bar:
        for hit in regions:
            ref = name2ref[hit.sacc]

            # sort the start/end (blast uses end<start for reverse strand)
            start, end = sorted((int(hit.sstart), int(hit.send)))

            result = _pileup.classic(bam, ref, start, end)

            if k_cor is not None:
                result.update(_pileup.experimental(bam, k_cor, kmer_length, fasta,
                                                   ref, start, end))

            if writer is None:
                fieldnames = ['sacc', 'start', 'end']
                fieldnames += sorted(list(result.keys()))
                writer = csv.DictWriter(outfile, fieldnames=fieldnames)
                writer.writeheader()

            result.update({
                'sacc': hit.sacc,
                'start': hit.sstart,
                'end': hit.send
            })
            writer.writerow(result)
            # bar.update(regionfile.tell())


@main.command()
@click.option("--readfile-type", "-t",
              type=click.Choice(['bam', 'fq']),
              help="Override filename based detection of readfile type.")
@click.option("--max-reads", "-m",
              type=click.IntRange(1), metavar="N",
              help="Only consider the first N reads.")
@click.option("--group-by", "-g",
              type=click.Choice(_scan.Flags.keys()), multiple=True, metavar="FLAG",
              help="Group output by BAM flag. May be specified multiple times. "
              "FLAG can be one of {}".format(list(_scan.Flags.keys())))
@click.option('--reference-fasta', '-f',
              type=click.File('rb'), metavar="FILE",
              help="Fasta file reads where mapped to. Required for boffset."
              " File may be bgzip'ed, but not gzip'ed.")
@click.option("--out-basehist", "-b",
              type=click.File("w", lazy=False), metavar="FILE",
              help="Compute histogram of base counts by position in read.")
@click.option('--boffset','-bo',
              type=click.IntRange(0,200), default=0, metavar="N", show_default=True,
              help="Include N bases prior to read start in base histogram. "
              "Requires reference fasta file.")
@click.option("--out-kmerhist", "-o",
              type=click.File("w", lazy=False), metavar="FILE",
              help="Compute histogram of kmers in reads")
@click.option("-k",
              type=click.IntRange(2,12), default=7, metavar="N", show_default=True,
              help="Length of kmers")
@click.option("--step", "-s",
              type=click.IntRange(1,100), default=7, metavar="N", show_default=True,
              help="Step between sampled kmers")
@click.option("--offset", "-O",
              type=click.IntRange(-100,100), default=0, metavar="N", show_default=True,
              help="Offset of first sampled kmer")
@click.option("--number", "-n",
              type=click.IntRange(1,100), default=8, metavar="N", show_default=True,
              help="Numer of sampled kmers")
@click.option('--out-mirrorhist', '-M',
              type=click.File('w', lazy=False), metavar="FILE",
              help="Compute histogram of mismatches against palindrome")
@click.option('--mirror-offset','-MO',
              type=click.IntRange(-100,100), default=4, metavar="N", show_default=True,
              help="Offset of palindrome center from read start")
@click.option('--mirror-length','-Ml',
              type=click.IntRange(1,50), default=10, metavar="N", show_default=True,
              help="Palindrome length is 2N+1")
@click.option('--out-isizehist', '-I',
              type=click.File('w', lazy=False), metavar="FILE",
              help="Compute histogram of insert sizes.")
@click.argument("readfile", nargs=-1, required=True)
def scan(readfile, readfile_type, out_basehist, boffset, out_kmerhist,
         k, number, step, offset, max_reads, reference_fasta,
         out_mirrorhist, mirror_offset, mirror_length,
         out_isizehist, group_by):
    """
    Gather read statistics
    """
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

    if readfile_type != 'bam' and reference_fasta:
        raise click.UsageError(
            "Reference fasta can only be used with mapped (bam/sam) reads"
        )

    update_every = 100000
    if max_reads and max_reads > 0 and max_reads/100 < update_every:
        update_every = int(max_reads+50/100)

    fasta = None

    if readfile_type == 'bam':
        infile = pysam.AlignmentFile(readfile[0])
        log.info("mapped = {}, unmapped = {}, total = {}"
                 "".format(
                     infile.mapped, infile.unmapped,
                     infile.mapped + infile.unmapped))
        length = infile.mapped + infile.unmapped
        if max_reads and max_reads > 0 and max_reads < length:
            length = max_reads

        def update_func(x):
            return x.update(update_every)

        if reference_fasta:
            fasta = pysam.FastaFile(reference_fasta.name)
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
        raise Exception("Internal error: Unknown readfile type?!")

    counters = []
    if out_basehist:
        counters.append(_scan.BaseHist(boffset))
    if out_kmerhist:
        counters.append(_scan.KmerHist(k, number, step, offset))
    if out_mirrorhist:
        counters.append(_scan.MirrorHist(mirror_offset, mirror_length))
    if out_isizehist:
        counters.append(_scan.IsizeHist())

    counters = _scan.ByFlag(counters,
                            [_scan.Flags[flag] for flag in group_by])

    with click.progressbar(length=length,
                           label="Scanning reads",
                           show_pos=True) as bar:
        with infile:
            nreads = _scan.scan_reads(
                infile, fasta, counters,
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

    if out_mirrorhist:
        out = csv.writer(out_mirrorhist)
        out.writerows(counters.get_rows(n))
        n = n + 1

    if out_isizehist:
        out = csv.writer(out_isizehist)
        out.writerows(counters.get_rows(n))
        n = n + 1


@main.command()
@click.option('--khist', '-k', type=click.File('r'), metavar="FILE",
              help="Match distribution of reads to statistics gathered "
                   "with `metacov scan`.")
@click.option('--klen', '-K', type=int, default=7,
              help="Length of k-mer")
@click.option('--proportions', '-p',
              help="Specify the proportions in which each genome should "
              "occur.")
@click.option("--out", "-o", "fwd_outfn", required=True, metavar="FILE",
              help="Name of output forward reads fastq file.")
@click.option("--out2", "-o2", "rev_outfn", required=True, metavar="FILE",
              help="Name of output reverse reads fastq file.")
@click.option("--num-reads", "-n", type=click.INT,
              help="Number of reads to sample")
@click.option("--rndSeed", "-rs", type=click.INT, default=1234,
              help="Random seed")
@click.option("--seqSys", "-ss",
              type=click.Choice(["GA1", "GA2", "HS10", "HS20", "HS25", "HSXn",
                                 "HSXt", "MinS", "MSv1", "MSv3", "NS50"]),
              default="HS20",
              help="Name of Illumina sequencing system")
@click.option("--len", "-l", "rlen", type=click.INT, default=100,
              help="Read length")
@click.option("--mflen", "-m", type=click.INT, default=450,
              help="Mean insert size")
@click.option("--sdev", "-s", type=click.INT, default=150,
              help="Stddev of insert size")
@click.argument("genomes", nargs=-1, required=True,
                metavar="FILE [FILE...]")
def simulate(khist, klen, proportions, fwd_outfn, rev_outfn,
             num_reads, rndseed, seqsys, rlen, mflen, sdev, genomes):
    """
    Pseudo-Randomly select reads from input file(s).

    Input read files may be in SAM, BAM or (gzipped) FastQ format
    (.sam, .bam, .fq, .fastq, .fq.gz or .fastq.gz).

    If a reverse read output filename is provided, input is assumed to be
    paired reads given either as sam/bam file(s) or pairs of FastQ files.
    """

    np.random.seed(rndseed)
    if proportions:
        proportions = [float(x) for x in proportions.split(",")]
        if len(proportions) != len(genomes):
            raise click.BadParameter(
                "List of proportions must be of equal length as list of"
                " genomes ({} != {})."
                "".format(len(proportions), len(genomes)),
                param_hint="-p")
    else:
        proportions = [1] * len(genomes)

    from metacov.simulate import generate_reads
    from contextlib import ExitStack

    def select_genome(proportions):
        proportions = np.array(proportions)
        proportions = proportions / np.sum(proportions)
        while True:
            vals = np.random.choice(len(proportions), size=255, p=proportions)
            for n in vals:
                yield n

    with ExitStack() as stack:
        outfq1 = stack.enter_context(FastQWriter(fwd_outfn))
        outfq2 = stack.enter_context(FastQWriter(rev_outfn))

        genome_sizes = []
        read_generators = []

        with click.progressbar(genomes, label="Launching ART") as bar:
            for genome in bar:
                size, readgen = stack.enter_context(
                    generate_reads(genome, seqsys, rlen, mflen, sdev, rndseed)
                )
                genome_sizes.append(size)
                read_generators.append(readgen)

        choose_genome = select_genome([x*y for x, y in zip(genome_sizes,
                                                           proportions)])

        if khist is not None:
            k_cor = _pileup.load_kmerhist(khist)
            for r in (0, 1):
                pmax = max(k_cor[r].values())
                k_cor[r] = {k: k_cor[r][k]/pmax for k in k_cor[r]}

        readno = 0
        used = 0

        bar = stack.enter_context(click.progressbar(
            length=num_reads, label="Generating Reads"))

        while True:
            used += 1
            genome_no = next(choose_genome)
            try:
                read = next(next(read_generators[genome_no]))
            except StopIteration:
                click.echo("Failed genome", genome_no, genomes[genome_no])
                click.echo("Produced {}/{} reads".format(readno, num_reads))
                break

            if khist:
                try:
                    k = (k_cor[0][read.read1.char_seq[:klen]] *
                         k_cor[1][read.read2.char_seq[:klen]])
                except KeyError:
                    # no match - usually means ambiguous base,
                    # skip this read (for simplicity reasons)
                    continue

                if np.random.sample() > k:
                    continue

            outfq1.write(read.read1)
            outfq2.write(read.read2)

            readno += 1
            bar.update(1)
            if readno >= num_reads:
                break

    click.echo("Used {} of {} reads (factor {})"
               "".format(used, readno, used/readno))
