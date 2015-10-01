#!/usr/bin/env python

import sys

num_regs = len(sys.argv) - 3

if sys.argv[1] in ['-h', '--help', '-help' '--h'] or num_regs < 1:
    print("Usage:- " + sys.argv[0] + " <matrix file> <bandwith> <positions (non negative)>")
    exit(0)

f = open(sys.argv[1])
f.readline()
vals = f.readline().split('\t')
bw = int(sys.argv[2])
for regs in range(0, num_regs):
    cReg = int(sys.argv[regs + 3])
    regAvg = 0.0
    for bp in range(cReg - bw, cReg + bw + 1):
        regAvg = regAvg + int(vals[bp + 1])
    regAvg = regAvg / ((bw * 2) + 1)
    print(str(cReg) + " +/- " + str(bw) + " avg = " + str(regAvg))
