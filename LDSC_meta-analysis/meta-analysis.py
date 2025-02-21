#!/usr/bin/env python

import gzip
import glob
import sys
import os
import tempfile
import subprocess
import math

from joblib import Parallel, delayed

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import matplotlib.font_manager



def dformat(xxl):
    x = xxl if type(xxl) is str else xxl[0]
    xx = xxl if type(xxl) is str else xxl[0]
    try:
        with open(x, 'rt') as f:
            if "LD Score Regression (LDSC)" not in f.read(): return
    except:
        return
    with open(x, 'rt') as f:
        lines = [ x.strip() for x in f.read().split("\n") ]

    lines = [ x if "/tmp" not in x else os.path.basename(xx).replace('.txt','') + "\t" + "\t".join(x.strip().split()[1:]) for x in lines ]
    
    
    if type(xxl) is list:
        def v(xx):
            print(xx, file = sys.stderr)
            with open(xx, 'rt') as f:
                llines = [ x.strip() for x in f.read().split("\n") ][-3:]
                llines = [ os.path.basename(xx).replace('.txt','') + "\t" + "\t".join(x.strip().split()[1:]) for x in llines if "/tmp" in x ]
                return llines
        r = Parallel(n_jobs = 64)(delayed(v)(xx) for xx in xxl[1:])
        for xx in r: lines += xx
    try:
        lidx = [ i for i, x in enumerate(lines) if x.startswith("Category") ][0]
    except:
        return
    with open(os.path.dirname(x) + "/formatted-h2/" + os.path.basename(x).replace('.txt','') + ".results", 'w') as o:
        o.write('\n'.join(lines[lidx:]))
    with open(os.path.dirname(x) + "/formatted-h2/" + os.path.basename(x).replace('.txt','') + ".log", 'w') as o:
        o.write('\n'.join(lines[:lidx]))

def run(x, d):
    with tempfile.NamedTemporaryFile('wt') as o:
        with open("formatted-h2.R", 'r') as f:
            o.write(f.read() % (x, d, ", ".join([ '"%s"' % os.path.basename(x).replace('.results','') for x in glob.glob(d + "/*.result*") ])))
        o.flush()
        with open("/dev/null", 'w') as dn:
            return subprocess.check_output("Rscript %s" % o.name, stderr = dn, shell = True).decode()

def kv(vx):
    r = {}
    lines = vx.strip().split("\n")
    for i in range(int(len(lines) / 2)):
        keys = [ x for x in lines[i * 2].strip().split() ]
        for ii, k in enumerate(keys):
            r[k] = float(lines[i * 2 + 1].strip().split()[ii + 1])
    return lines[1].strip().split()[0], r

def trun(x, d):
    try:
        return kv(run(x, d))
    except:
        return None

def percentages(inputf):
    common = glob.glob("/data/zusers/sheddn/LDR/common_snps/*common*bed")
    def read_c(cc):
        with open(cc, 'r') as f:
            return { x.strip().split()[-1] for x in f }
    r = Parallel(n_jobs = 64)(delayed(read_c)(x) for x in common)
    all_common = set()
    for x in r:
        all_common = all_common.union(x)
    def totals(iif):
        count = 0
        with gzip.open(iif, 'rt') as ff:
            hmap = { i: x.strip() for i, x in enumerate(ff.readline().split("\t")) }
            counts = { x: 0 for _, x in hmap.items() }
            for line in ff:
                if line.strip().split()[2] not in all_common: continue
                count += 1
                for i, x in enumerate(line.strip().split()):
                    counts[hmap[i]] += 1 if x == "1" else 0
        return count, counts
    results = Parallel(n_jobs = 64)(delayed(totals)(x) for x in inputf)
    gtotal = sum([ x[0] for x in results ])
    return { k.strip().replace('[', "").replace(']', "").replace(',', ""): sum([ x[1][k] if len(x) >= 2 and k in x[1] else 0 for x in results ]) / float(gtotal) for k, _ in results[0][1].items() }

def main(argc, argv):

    if argc < 3:
        print("usage: meta-analysis.py ldsc-results-directory *.annot.gz", file = sys.stderr)
        return 1

    os.system("mkdir -p %s" % (argv[1] + "/formatted-h2"))
    for x in glob.glob(argv[1] + "/*"):
        # try:
        #     dformat(x)
        # except:
        #     print(x)
        #     continue
        dformat(x)
    print("running meta analysis for %d traits" % len(glob.glob(argv[1] + "/formatted-h2/*results")), file = sys.stderr)
    values = Parallel(n_jobs = 24)(delayed(trun)(x, argv[1] + "/formatted-h2") for x in range(2, 500))

    cmap = { x[0]: x[1] for x in values if x is not None }
    p = percentages(argv[2:])
    print("partition\tpercent SNPs\tenrichment\terror\tp\ttau\ttau error\ttau p")
    for k, v in cmap.items():
        # print("%s\t%f\t%f\t%f\t%e\t%f\t%f\t%e" % (k.replace("L2_0", ""), p[k.split("L2")[0]], v["Enr"], v["Enr_se"], v["Enr_P"], v["tau"], v["tau_se"], v["tau_P"]))
        print("%s\t%f\t%f\t%f\t%e\t%f\t%f\t%e" % (k.replace("L2_0", ""), p[k.split("L2_0")[0]], v["Enr"], v["Enr_se"], v["Enr_P"], v["tau"], v["tau_se"], v["tau_P"]))
    return 0

if __name__ == "__main__":
    sys.exit(main(len(sys.argv), sys.argv))
