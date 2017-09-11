import click
import logging
import metacov.blast
import pysam
import csv
import os
from io import open

from metacov import pileup

logging.basicConfig(
    level=logging.INFO,
    format="[%(relativeCreated)6.1f %(funcName)s]  %(message)s",
    datefmt="%I:%M:%S"
)

log = logging.getLogger(__name__)


@click.command(short_help="""
BAMFILE      Input BAM file. Must be sorted and indexed
REGEIONFILE  Input Region file. Can be one of: BLAST 7
""")
@click.argument('bamfile', type=click.File('rb'))
@click.argument('regionfile_name')#, type=click.File('r'))
@click.argument('coveragefile', type=click.File('w'))
def main(bamfile, regionfile_name, coveragefile):
    """
    Metacov -- Calculate abundance estimates over metagenome regions
    
    \b
    Arguments:
      BAMFILE      Input BAM file. Must be sorted and indexed
      REGEIONFILE  Input Region file. Can be one of: BLAST 7
    """
    regionfile = open(regionfile_name, "r")
    log.info("Gathering coverages for regions in \"{}\"".format(regionfile.name))
    blastreader = blast.reader(regionfile)        


    log.info("Processing BAM file \"{}\"".format(bamfile.name))
    bam = pysam.AlignmentFile(bamfile.name)

    try:
        log.info(("Number of reads:\n"
                  "  total:    {total}\n"
                  "  mapped:   {mapped} ({mappct}%)\n"
                  "  unmapped: {unmapped}\n"
        ).format(
            mapped=bam.mapped, unmapped=bam.unmapped,
            total=bam.mapped+bam.unmapped,
            mappct=bam.mapped / (bam.mapped+bam.unmapped) * 100
        ))
    except:
        log.error("BAM file not indexed!?")

    name2ref = { word.split()[0]:word
                 for word in bam.references }

    writer = None

    regionfile_size = os.path.getsize(regionfile.name)
    with click.progressbar(
#            length=os.path.getsize(regionfile.name),
            length=0,
            label="Calculating coverages") as bar:
        for hit in blastreader:
            # translate hit name, in case the bam file contains
            # multi-word references
            ref = name2ref[hit.sacc]

            # sort the start/end (blast uses end<start for reverse strand)
            start, end = sorted((hit.sstart, hit.send))
            length = end-start

            classic = pileup.classic(bam, ref, start, end)

            if writer is None:
                fieldnames = [ 'sacc', 'start', 'end' ] + sorted(list(classic.keys()))
                writer = csv.DictWriter(coveragefile, fieldnames=fieldnames)
                writer.writeheader()

            classic.update({
                'sacc': hit.sacc,
                'start':hit.sstart,
                'end': hit.send
            })
            writer.writerow(classic)
            #bar.update(regionfile.tell())



        
                      

    
    
    
