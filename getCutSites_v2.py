#!/usr/bin/env python

import pysam
import argparse
import sys
import time
import collections as c
import math
import getpass
import socket
import re

if socket.gethostname() == 'GBW-S-SEQ03':
    print("Working on SEQ-SRV-03 - fixing path...")
    for n, i in enumerate(sys.path):
        if re.match(".*pandas.*", i) is not None:
            del(sys.path[n])

import pandas as pd


__author__ = 'kdavie'

""" Takes a sam or bam file as an input and outputs a 6 column bed file
with the locations of each read start point and its strand, col 5 = .
Also able to create a matrix of either read depth or fragment start sites
based off of a bed input. When the matrix option is used, a matrix will
be created based on the defined parameters, this can easily be used to
create heatmaps and aggregation plots.
"""

# Help and arguments parser
parser = argparse.ArgumentParser(description='Takes a sam or bam file as an input and outputs a 6 column bed file '
                                             'with the locations of each read start point and its strand, col 5 = . '
                                             'Also able to create a matrix of either read depth or fragment start sites '
                                             'based off of a bed input. When the matrix option is used, a matrix will '
                                             'be created based on the defined parameters, this can easily be used to '
                                             'create heatmaps and aggregation plots.')
parser.add_argument('inputFile', metavar='input', type=str, help='Input data file - Either SAM, BAM or Tabix formatted '
                                                                 'bed file (autodetected)', default=sys.stdin)
parser.add_argument('bed', metavar='bed', type=str, help='Bed file for creating a matrix')
parser.add_argument('size', metavar='size', type=int, help='How wide should we make the matrix from each peak?',
                    default=0)
parser.add_argument('outputFile', metavar='output', type=str, help='Output file - either bed file or matrix',
                    default=sys.stdout)
parser.add_argument('--quiet', action='store_true', help='Outputs no information || false', default=False)
parser.add_argument('-p', type=float, help='Percentage increment that time and progress is reported for || 10',
                    default=10)
parser.add_argument('--atac', action='store_true', help='Processes sites with a 9bp difference to account for binding'
                                                        ' vs cutting in ATAC-seq || false', default=False)
parser.add_argument('-pos', metavar='bp', type=int, help='Override default ATAC-seq positive shift || 4',
                    default=0)
parser.add_argument('-neg', metavar='bp', type=int, help='Override default ATAC-seq Negative shift || -5',
                    default=0)
parser.add_argument('--debug', action='store_true', help='Print debug info, confirm every read', default=False)
parser.add_argument('-smooth', metavar='bp', type=int, help='Smoothening to apply in bp, only for cutsites', default=1)
parser.add_argument('-binbed', action='store_true', help='Process bed as binary?')
parser.add_argument('--sumsOnly', action='store_true', help='Only provide the sums of the matrix')
parser.add_argument('--rpm', action='store_true', help='Normalize the matrix to reads per million')
parser.add_argument('--regNorm', action='store_true', help='Normalize the matrix to number of regions')
parser.add_argument('-ylim', type=float, help='Set the upper y limit (lower will be 0)')
args = parser.parse_args()

if args.debug:
    print
    print
    print("Debug info:")
    print
    print(" Input file:             " + args.inputFile)
    print(" Output file:            " + args.outputFile)
    print(" Bed file used           " + str(args.bed))
    print(" Quiet?:                 " + str(args.quiet))
    print(" Reporting increment:    " + str(args.p) + "%")
    print(" ATAC-Seq?               " + str(args.atac))
    print(" Manual Positive shift   " + str(args.pos))
    print(" Manual Negative shift   " + str(args.neg))
    print(" Smoothing window        " + str(args.smooth))
    print(" Bed vs Bed as binary?   " + str(args.binbed))
    print
    print

# Try and load the data file - Can work with both B/SAM and Tabix formats
try:
    datafile = pysam.Samfile(args.inputFile, "rb")
    isTabix = False
except ValueError:
    try:
        datafile = pysam.Tabixfile(args.inputFile, "r")
        isTabix = True
        if not args.quiet:
            print("Input is bed format")
    except ValueError:
        print("Input file is of an unsupported type")
        exit(1)
except OSError:
    print("ERROR: File not found!")
    exit(1)

try:
    b = open(args.bed, 'r')
except TypeError:
    print("Cannot open bed file")
    exit(1)
except OSError:
    print("ERROR: File not found!")
    exit(1)

# Assign some initial variables
count = 0
outBed = []
reverse = 0
forward = 0
percent = 0
bufferPos = 0
bufferNeg = 0
bedLen = 0
zerr = 0
mSums = ["Sums"]
if args.rpm and isTabix:
    print("Cannot RPM normalize a bed file - will continue without normalisation")
    args.rpm = False
if args.rpm:
    libSize = datafile.mapped

# Deal with ATAC overrides
if args.atac:
    if args.bufferPos > 0:
        bufferPos = args.pos
    else:
        bufferPos = 4
    if args.bufferNeg > 0:
        bufferNeg = args.neg
    else:
        bufferNeg = -5

# Function for making a matrix
def makematrix(n, region):
    global count, percent, totalLen, end_time, start_time, readName, forward, reverse, outBedNew, outBed, bufferPos, \
        bufferNeg, f, bedLen, isTabix, zerr

    # Split the input bed
    if args.debug:
        print("Working on - " + region)
    chrom = region.split('\t')[0]
    start = int(region.split('\t')[1])
    end = int(region.split('\t')[2])
    reg = str(line[0]) + ':' + str(line[1]) + '-' + str(line[2])
    count += 1
    cnt = c.defaultdict(lambda: 0)
    cntneg = c.defaultdict(lambda: 0)
    totalreads = 0

    # If a size is specified change the start and end sites to be Half way between +/- size, if not
    oldstart = start
    oldend = end

    start = int((oldstart + (math.floor((oldend - oldstart) / 2))) - args.size)
    end = int((oldstart + (math.floor((oldend - oldstart) / 2))) + args.size)
    matrixline = []

    # Only process reads which are within the region boundaries
        # Report Percentage based on current region vs total
    try:
        if count % round((bedLen / (100.0 / args.p))) == 0 and not args.quiet:
            percent += args.p
            end_time = time.time()
            print("Processed " + str(percent) + "% (" + str(count) + " of " + str(bedLen) + " regions) in %.2f"  %  (end_time - start_time) + \
            " seconds.")
            start_time = time.time()
    except ZeroDivisionError:
        if zerr == 0:
            print("Too few regions, will not report progress")
            zerr = 1

    # Try and extract strand information, if non existent or not '+'/'-', assume positive strand
    try:
        strand = region.split('\t')[5].rstrip('\n')
        if strand != '+' or strand != '-':
            strand = '+'
    except IndexError:
        strand = '+'

    startpileup = time.time()
    # Support for Tabix indexed bed files
    if isTabix:
        try:
            for reg in datafile.fetch(chrom, start, end):
                sreg = reg.split('\t')
                if sreg[0] == chrom and (start <= int(sreg[1]) <= end or
                                         start <= int(sreg[2]) <= end):
                    for i in range(int(sreg[1]), int(sreg[2]), 1):
                        if args.binbed:
                            cnt[i] += 1
                        try:
                            cnt[i] += int(float(sreg[4]))
                        except (IndexError, ValueError):
                            cnt[i] += 1
        except ValueError:
            if args.debug:
                print("Region " + chrom + " " + str(start) + " " + str(end) + " does not exist in this file")
        if args.debug:
            print("Finished checking bed in " + str(time.time() - startpileup) + " seconds")

        if strand == '+':
            for base in range(start, end + 1, 1):
                if base < start or base > end:
                    continue
                totalreads += int(cnt[base])
                matrixline.append(float(cnt[base]).rstrip('\n'))

        if strand == '-':
            for base in range(end + 1, start, 1):
                if base < start or base > end:
                    continue
                totalreads += int(cnt[base])
                matrixline.append(float(cnt[base]).rstrip('\n'))
            matrixline = matrixline[::-1]
    # Support for SAM and BAM files
    if not isTabix:
        cols = 0
        nextcol = start
        if start < 0:
            start = 0
        for pileupColumn in datafile.pileup(chrom, start, end):
            # if args.debug:
                # print(start, pileupColumn.pos, end)
            if int(start) > int(pileupColumn.pos) or int(end) < int(pileupColumn.pos):
                continue
            while nextcol < pileupColumn.pos:
                matrixline.append(0)
                nextcol += 1
                cols += 1
            matrixline.append(float(pileupColumn.n))
            if args.debug:
                print("Appending " + str(float(pileupColumn.n)))
            totalreads += pileupColumn.n
            nextcol += 1
            cols += 1
        if args.debug:
            print("Total number of pileup columns: " + str(cols))
            print("Command sent: reads.pileup(" + str(chrom) + ", " + str(start) + ", " + str(end) + ")")
            print("cols - (end - start) = " + str(cols - (end - start)))
            print(str(nextcol) + ", " + str(end) + ", " + str(start) + ", " + str(cols))
        if nextcol <= end:
            for i in range(nextcol, end + 1, 1):
                matrixline.append(0)
        if strand == '-':
            matrixline = matrixline[::-1]
        return(matrixline)


if not args.quiet:
    print("\nBeginning to process file\n")

# Feedback and run the main loops

index = []
bedLines = []
for line in b:
    line = line.rstrip('\n').split('\t')
    [str(x) for x in line]
    lineName = line[0] + ':' + line[1] + '-' + line[2]
    if len(line) >= 3:
        if lineName not in index:
            bedLen += 1
            bedLines.append("\t".join(line))
            index.append(lineName)
        else:
            if not args.quiet:
                print("Duplicate entry found - " + "\t".join(line))
    elif len(line) >= 1:
            if not args.quiet:
                print("Malformed line - " + "\t".join(line))

cols = [x for x in range(-(args.size), args.size + 1)]
# matrix = pd.DataFrame(index=index, columns=cols)

mlines = []
start_time = time.time()
for n, r in enumerate(bedLines):
    newLine = makematrix(n, r)
    if newLine and len(newLine) == ((args.size * 2) + 1):
        mlines.append(newLine)

if not args.quiet:
    print("Generating matrix...")
matrix = pd.DataFrame(mlines, index=index, columns=cols)

if args.rpm:
    if not args.quiet:
        print("Normalising matrix...")
    matrix = (matrix.fillna(0) / libSize) * 1000000

if args.regNorm:
    matrix = matrix / int(bedLen)

if args.smooth > 1:
    if not args.quiet:
        print("Smoothing matrix...")
    matrix = pd.rolling_mean(matrix.T, args.smooth).T
    if not args.quiet:
        print("Selecting columns")
    smoothMatrix = pd.DataFrame(index=index)
    for i in matrix.columns:
        if args.debug:
            print("On column" + str(i))
        if int(i) % args.smooth == 0 or int(i) == 0:
            if int(i) == args.size or int(i) == -(args.size):
                continue
            smoothMatrix[i] = matrix[i]
    matrix = smoothMatrix

if not args.quiet:
    print("Calculating sums...")

matrixSums = pd.DataFrame(matrix.sum()).T


if not args.quiet:
    print("Generating plot...")


import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

if args.smooth > 1:
    matrix.sum().plot(xlim=(-(int(args.size)) + 10, (int(args.size)) - 10), color="black")
else:
    matrix.sum().plot(xlim=(-(int(args.size)), (int(args.size))), color="black")

if args.rpm:
    plt.ylabel('Aggregate reads per (mapped) million')
else:
    plt.ylabel('Aggregate reads')

plt.xlabel('Position')

if args.ylim is not None:
    plt.ylim(0, args.ylim)
plt.axvline(x=0, ls='--', lw=1)
plt.savefig('.'.join([args.outputFile, 'pdf']), bbox_inches='tight')

if args.rpm:
    matrixSums.index = ['Aggregate Reads per (mapped) million']
else:
    matrixSums.index = ['Aggregate Reads']
if args.sumsOnly:
    matrixSums.to_csv(args.outputFile, sep='\t', na_rep='0')
else:
    if not args.quiet:
        print("Writing full matrix. This can take some time and a lot of space.")
    matrix.to_csv(args.outputFile, sep='\t', na_rep='0')
