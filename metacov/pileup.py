import pysam
import numpy as np

def classic(bam, ref, start, end):
    length = end - start
    columns = np.zeros(length)

    for column in bam.pileup(ref, start, end):
        if column.pos < start or column.pos >= end:
            continue
        columns[column.pos - start] += column.n

    return {
        'min': np.amin(columns),
        'max': np.amax(columns),
        'med': np.median(columns),
        'std': np.std(columns),
        'avg': np.mean(columns),
        'q23': np.mean(sorted(columns)[length//4:length-length//4]),
        'sum': np.sum(columns)
    }

