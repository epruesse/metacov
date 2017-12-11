import os
import subprocess as sp
from contextlib import contextmanager, ExitStack
from tempfile import TemporaryDirectory

from metacov.pyfq import FastQFilePair


@contextmanager
def generate_reads(genome_fn, model, readlen, fraglen, fragsd,
                   seed=1234, nf=5):
    if sp.call(["command", "-V", "art_illumina"], stdout=sp.DEVNULL, shell=True) != 0:
        raise Exception("Command not found: art_illumina")
    with ExitStack() as stack:
        tmpdir = stack.enter_context(TemporaryDirectory())
        outfq1_fn = os.path.join(tmpdir, "out1.fq")
        outfq2_fn = os.path.join(tmpdir, "out2.fq")
        infa_fn = os.path.join(tmpdir, "in.fasta")
        if genome_fn.endswith(".gz"):
            sp.call("gunzip -c {} > {}".format(genome_fn, infa_fn), shell=True)
        else:
            os.link(genome_fn, infa_fn)

        genome_len = 0
        with open(infa_fn) as infa:
            for line in infa:
                if line.startswith(">"):
                    continue
                genome_len += len(line) - line.count(' ')

        os.mkfifo(outfq1_fn)
        os.mkfifo(outfq2_fn)

        stack.enter_context(sp.Popen([
                    "art_illumina",
                    "--seqSys", model,
                    "-nf", str(nf),
                    "--noALN",
                    "--in", infa_fn,
                    "--len", str(readlen),
                    "--fcov", str(1000000000),
                    "--mflen", str(fraglen),
                    "--sdev", str(fragsd),
                    "--out", os.path.join(tmpdir, "out"),
                    "--paired",
                    "--rndSeed", str(seed)
        ], shell=False, stdout=sp.DEVNULL, stderr=sp.DEVNULL))

        fqfp = stack.enter_context(FastQFilePair(outfq1_fn, outfq2_fn))
        yield (genome_len, fqfp)
