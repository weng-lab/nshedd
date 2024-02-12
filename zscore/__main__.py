#!/usr/bin/env python

import sys

sys.dont_write_bytecode = True

import gzip
import statistics
import pyBigWig
import pysam
import math
import numpy

from functools import reduce
from joblib import Parallel, delayed

def batch(iterable, n=1):
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx:min(ndx + n, l)]

def flatten(x):
    r = []
    for xx in x:
        r += xx
    return r
        
def read(bigwig, regions):
    with pyBigWig.open(bigwig) as bw:
        return [ statistics.mean(
            bw.values(x.split()[0], int(x.split()[1]), int(x.strip().split()[2]))
        ) for x in f ]

def tmean(bw, *args):
    try:
        return statistics.mean(bw.values(*args))
    except:
        return 0.0

def contains(read, start, end):
    position = read.reference_start if not read.is_reverse else read.reference_end
    return position is not None and position > start and position < end

def bamcount(chromosome, start, end, bam, total, i = None):
    if i is not None: print(i, file = sys.stderr)
    return len([ x for x in bam.fetch(chromosome, start, end) if contains(x, start, end) ])

def bamzscoress(regions, bam, total):
    b = pysam.AlignmentFile(bam, 'rb')
    r = [ bamcount(x.split()[0], int(x.split()[1]), int(x.strip().split()[2]), b, total, "%d / %d" % (i, len(regions)) if i % 1000 == 0 else None ) for i, x in enumerate(regions) ]
    b.close()
    return r

def bamzscores(bed, bam, j = 1):
    with (gzip.open if bed.endswith(".gz") else open)(bed, 'rt') as f:
        lines = [ x.strip() for x in f ]
    total = sum([ int(l.strip().split('\t')[2]) for l in pysam.idxstats(bam) if len(l.strip().split('\t')) >= 3])
    bsetlen = int(math.ceil(len(lines) / float(j)))
    breakpoints = [ bsetlen * x for x in range(j + 1) ]
    jsets = [ lines[breakpoints[i] : breakpoints[i + 1]] for i in range(len(breakpoints) - 1) ]
    jresults = Parallel(n_jobs = j)(delayed(bamzscoress)(x, bam, total) for x in jsets)
    allresults = []
    for x in jresults:
        allresults += x
    return ztpm(allresults), sum(allresults)

def ztpm(signal):
    treads = sum(signal)
    zeros = len([ x for x in signal if x == 0 ])
    r = [ x * 100000 / float(treads + zeros) for x in signal ]
    m = numpy.mean([ math.log(x) for x in r if x > 0 ])
    s = numpy.std([ math.log(x) for x in r if x > 0 ])
    return [ (math.log(x) - m) / s if x > 0 else -10 for x in r ]

def zt(signal):
    signalmean = statistics.mean([ math.log(x) for x in signal if x > 0.0 ])
    signalstd = statistics.stdev([ math.log(x) for x in signal if x > 0.0 ])
    return [ (math.log(x) - signalmean) / signalstd if x > 0.0 else -10.0 for x in signal ]

def tmeans(bigwig, group):
    bw = pyBigWig.open(bigwig)
    if bw is None:
        raise Exception("Error opening %s: no such file or directory." % bigwig)
    return [ tmean(bw, x.split()[0], int(x.split()[1]), int(x.strip().split()[2])) for x in group ]

def zscores(bed, bigwig, j = 1):
    with (gzip.open if bed.endswith(".gz") else open)(bed, 'rt') as f:
        l = [ x for x in f ]
    signal = flatten(Parallel(n_jobs = 48)(delayed(tmeans)(bigwig, x) for x in batch(l, 48)))
    return zt(signal), None

def main(argc, argv):

    if argc < 3:
        print("usage: zscore signal.{bigwig|bam} peaks.bed [threshold] [j] [number-of-reads-in-peaks.txt]", file = sys.stderr)
        return 1

    z, t = (zscores if not argv[1].endswith(".bam") else bamzscores)(argv[2], argv[1], int(argv[4]) if argc >= 5 else 1)
    if argc >= 6:
        with open(argv[-1], 'r') as f:
            sortedscores = sorted([ float(x.strip().split()[-1]) for x in f ])
        sortedindexes = sorted(range(len(z)), key = lambda i: z[i])
        sortedindexes = { v: k for k, v in enumerate(sortedindexes) }
    else:
        sortedscores = z
        sortedindexes = range(len(z))
    with (gzip.open if argv[2].endswith(".gz") else open)(argv[2], 'rt') as f:
        for idx, line in enumerate(f):
            if argv[3]=="NA":
                print("%s\t%f" % (line.strip(), sortedscores[sortedindexes[idx]]))
            else:
                if sortedscores[sortedindexes[idx]] > float(argv[3]):
                    print("%s\t%f" % (line.strip(), sortedscores[sortedindexes[idx]]))
    if argv[1].endswith(".bam") and argc >= 5:
        with open(argv[-1], 'w') as o:
            o.write(str(t) + '\n')

    return 0

if __name__ == "__main__":
    sys.exit(main(len(sys.argv), sys.argv))
