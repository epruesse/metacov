import numpy as np

import pandas as pd

import scipy.stats as sps

import sys

def classic(bam, ref, start, end):
    length = end - start
    columns = np.zeros(length)

    for column in bam.pileup(ref, start, end):
        if column.pos < start or column.pos >= end:
            continue
        columns[column.pos - start] += column.n

    return {
        'min': int(np.amin(columns)),
        'max': int(np.amax(columns)),
        'med': int(np.median(columns)),
        'std': round(np.std(columns), 2),
        'avg': round(np.mean(columns), 2),
        'q23': round(np.mean(sorted(columns)[length//4:length-length//4]), 2),
        'sum': int(np.sum(columns))
    }


def load_kmerhist(f):
    df = pd.read_csv(f)
    df.drop(df[(df.Mapped == "Unmapped") | (df.kmer == "N"*7)].index,
            inplace=True)
    df.set_index('kmer', inplace=True)
    cor = df[df.columns[0]] / df[df.columns[1:]].mean(axis=1)
    return [cor[df.R == r].to_dict() for r in ("R1", "R2")]


def experimental(bam, k_cor, fasta, ref, start, end):
        length = end-start
        if length == 0:
            raise Exception("Length must be > 0")
        cov = np.zeros(length)
        cov2 = np.zeros(length)
        wnf = 0.0
        cov_cor = np.zeros(length)
        cor = np.zeros(length)
        starts = np.zeros(length)
        nreads = 0
        i = 0
        secondary = 0
        improper = 0
        x = {}

        insert = 450
        sd = 150
        #istart = insert - 2 * sd
        #iend = insert + 2 * sd
        istart = 0
        iend = insert * 2
        norm = sps.norm(insert, sd).pdf(range(istart, iend + 1)).ravel()

        if fasta:
            region = fasta.fetch(ref, start, end).upper()
            gc = region.count("G") + region.count("C")
            gc = gc / (gc + region.count("A") + region.count("T"))
            if k_cor:
                cor_fwd = np.zeros(length)
                cor_rev = np.zeros(length)
                n = 7
                for i in range(length - n):
                    try:
                        cor_fwd[i] = k_cor[0][region[i:i+n]]
                    except KeyError:
                        cor_fwd[i] = 0
                for i in range(-length + n - 1, 0):
                    try:
                        cor_rev[length + i] = k_cor[1][region[i:i-n:-1]]
                    except KeyError:
                        cor_rev[length + i] = 0

                cor_revsum = np.zeros(length)
                for i in range(length):
                    l = min(length - i, iend - istart)
                    cor_revsum[i] = np.inner(norm[0:l], cor_rev[i:i+l])
                ecor = np.inner(cor_fwd, cor_revsum) / length
        else:
            region = None
            gc = "na"
            ecor = "na"

        for read in bam.fetch(ref, start, end):
            # skip but count secondary reads
            if read.is_secondary:
                secondary += 1
                continue

            # skip improper reads
            if not read.is_proper_pair:
                improper += 1
                continue

            if read.query_name in x:
                readb = x[read.query_name]
                s = min(read.reference_start, readb.reference_start) - start
                e = max(read.reference_start, readb.reference_start) - start
                cov2[s-1:e+1] += 1
                try:
                    wnf += 1/(
                        k_cor[1 if read.is_reverse else 0][read.query_alignment_sequence[0:7]]
                        *
                        k_cor[1 if readb.is_reverse else 0][readb.query_alignment_sequence[0:7]]
                    )
                except KeyError:
                    wnf += 1
                except ZeroDivisionError:
                    wnf += 1
                del x[read.query_name]
            else:
                x[read.query_name] = read

            # if not read.reference_end: continue

            readno = 0 if read.is_read1 else 1
            kmer = read.query_alignment_sequence[0:7]
            try:
                rcor = k_cor[readno][kmer]
            except:
                rcor = 1
            if rcor == 0:
                print("RCOR is ZERO: {} {}".format(["R1", "R2"][readno], kmer))
                rcor = 1

            if read.is_reverse:
                rend = read.reference_start - start
                rstart = rend - read.reference_length
            else:
                rstart = read.reference_start - start
                rend = rstart + read.reference_length

            for i in range(max(0, rstart), min(length, rend)):
                cov[i] += 1
                cov_cor[i] += (1/rcor)

            if rstart >= 0 and rstart < length:
                cor[rstart] = 1/rcor
                starts[rstart] = 1
                nreads += 1
                # rpkm += 1
                # rpkm_cor += 1/rcor

        nz = sum([1 for n in starts if n == 0])
        nz_e = length * (1 - 1/length)**nreads
        nzef = nz / nz_e

        # rpkm /= length / 1000 * tot_reads / 1000000
        # rpkm_cor /= length / 1000 * tot_reads / 1000000

        allreads = secondary+nreads+improper

        return {
            'cov':   np.mean(cov),
            'covc':  np.mean(cov_cor),
            'den':   round(np.mean(starts), 3),
            'denc':  round(np.mean(cor), 3),
            'cov2':  round(np.mean(cov2)),
            'cf':    round(sum(cor)/sum(starts), 3),
            'ambig': round(secondary/(allreads), 3) if allreads > 0 else 0,
            'improper': round(improper/(allreads), 3) if allreads > 0 else 0,
            'nzef':  round(nzef, 3),
            'gc': round(gc, 3),
            'ecor' : round(ecor, 3),
            'wnf': round(wnf/length, 3),
            'cov3': round(200 * (wnf/ecor) / nzef / length,3)
        }
