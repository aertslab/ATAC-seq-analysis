#!/usr/bin/env python

import pysam
import argparse
import sys
import time
import collections as c
import math

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
parser.add_argument('outputFile', metavar='output', type=str, help='Output file - either bed file or matrix',
                    default=sys.stdout)
parser.add_argument('--addChr', action='store_true', help='Prepends chr to all lines || false', default=False)
parser.add_argument('--quiet', action='store_true', help='Outputs no information || false', default=False)
parser.add_argument('-p', type=float, help='Percentage increment that time and progress is reported for || 10',
                    default=10)
parser.add_argument('--atac', action='store_true', help='Processes sites with a 9bp difference to account for binding'
                                                        ' vs cutting in ATAC-seq || false', default=False)
parser.add_argument('-pos', metavar='bp', type=int, help='Override default ATAC-seq positive shift || 4',
                    default=0)
parser.add_argument('-neg', metavar='bp', type=int, help='Override default ATAC-seq Negative shift || -5',
                    default=0)
parser.add_argument('--matrix', action='store_true', help='Instead of making a bed file, a matrix is created || false',
                    default=False)
parser.add_argument('--pileup', action='store_true', help='Count reads instead of cutSites || False', default=False)
parser.add_argument('-size', metavar='bp', type=int, help='How wide should we make the matrix from each peak?',
                    default=0)
parser.add_argument('-bed', metavar='file', type=str, help='Bed file for creating a matrix')
parser.add_argument('--debug', action='store_true', help='Print debug info, confirm every read', default=False)
parser.add_argument('--centipede', action='store_true', help='Required for a correct matrix for centipede (ALPHA)',
                    default=False)
parser.add_argument('-smooth', metavar='bp', type=int, help='Smoothening to apply in bp, only for cutsites', default=1)
parser.add_argument('-binbed', action='store_true', help='Process bed as binary?')
parser.add_argument('--sumsOnly', action='store_true', help='Only provide the sums of the matrix')
args = parser.parse_args()

if args.debug:
    print
    print
    print("Debug info:")
    print
    print(" Input file:             " + args.inputFile)
    print(" Output file:            " + args.outputFile)
    print(" Bed file used           " + str(args.bed))
    print(" Adding chr?:            " + str(args.addChr))
    print(" Quiet?:                 " + str(args.quiet))
    print(" Reporting increment:    " + str(args.p) + "%")
    print(" ATAC-Seq?               " + str(args.atac))
    print(" Manual Positive shift   " + str(args.pos))
    print(" Manual Negative shift   " + str(args.neg))
    print(" Outputting matrix?      " + str(args.matrix))
    print(" Doing pileup?           " + str(args.pileup))
    print(" Centipede compatible?   " + str(args.centipede))
    print(" Smoothing window        " + str(args.smooth))
    print(" Bed vs Bed as binary?   " + str(args.binbed))
    print
    print

if args.matrix and not args.bed:
    print("A bed file is required to create a matrix")
    exit()

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

if not args.matrix and isTabix:
    print("Cannot make a bed of a bed")
    exit(1)

f = open(args.outputFile, 'w')
if args.bed:
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
pairedName = set()
start_time = time.time()
percent = 0
bufferPos = 0
bufferNeg = 0
bedLen = 0
zerr = 0
mSums = ["Sums"]

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


def getsite(i):
    global count, percent, totalLen, end_time, start_time, readName, forward, reverse, outBedNew, outBed, bufferPos, \
        bufferNeg, f
    count += 1

    # Tid is Template ID or chromosome name
    if i.tid > -1:
        try:
            if count % round((totalLen / (100.0 / args.p))) == 0 and not args.quiet:
                percent += args.p
                end_time = time.time()
                print("Processed " + str(percent) + "% (" + str(count) + " of " + str(totalLen) + " reads) in %.2f"  %  (end_time - start_time) + \
                " seconds.")
                start_time = time.time()
        except ZeroDivisionError:
            if zerr == 0:
                print("Too few regions, will not report progress")
                zerr = 1

        # Save read names for later
        if str(i.qname) not in pairedName:
            pairedName.add(i.qname)
            readName = str(i.qname) + "_1"
        else:
            readName = str(i.qname) + "_2"
            pairedName.remove(i.qname)

        # Deal with forward and reverse reads separately
        if not i.is_reverse:
            forward += 1

            # If we are pre-pending chr, don't prepend to the following chromosome names
            if args.addChr:
                if str(datafile.getrname(i.tid)) == 'dmel_mitochondrion_genome':
                    outBedNew = str(datafile.getrname(i.tid)) + '\t' + str(i.pos + bufferPos) + '\t' + str(
                        i.pos + bufferPos) + '\t' + readName + '\t' + "." + '\t' + "+"
                else:
                    outBedNew = 'chr' + str(datafile.getrname(i.tid)) + '\t' + str(i.pos + bufferPos) + '\t' + str(
                        i.pos + bufferPos) + '\t' + readName + '\t' + "." + '\t' + "+"
            else:
                outBedNew = str(datafile.getrname(i.tid)) + '\t' + str(i.pos + bufferPos) + '\t' + str(
                    i.pos + bufferPos) + '\t' + readName + '\t' + "." + '\t' + "+"

        # If the region is on the reverse strand, we get the last base of the read instead of the first
        # Inferred_length uses the CIGAR score to predict the length of the read.
        if i.is_reverse:
            reverse += 1
            if args.addChr:
                if str(datafile.getrname(i.tid)) == 'dmel_mitochondrion_genome':
                    outBedNew = str(datafile.getrname(i.tid)) + '\t' + str(
                        i.pos + i.inferred_length + bufferNeg) + '\t' + str(
                        i.pos + i.inferred_length + bufferNeg) + '\t' + readName + '\t' + "." + '\t' + "-"
                else:
                    outBedNew = 'chr' + str(datafile.getrname(i.tid)) + '\t' + str(
                        i.pos + i.inferred_length + bufferNeg) + '\t' + str(
                        i.pos + i.inferred_length + bufferNeg) + '\t' + readName + '\t' + "." + '\t' + "-"
            else:
                outBedNew = str(datafile.getrname(i.tid)) + '\t' + str(
                    (i.pos + i.inferred_length + bufferNeg)) + '\t' + str(
                    i.pos + i.inferred_length + bufferNeg) + '\t' + readName + '\t' + "." + '\t' + "-"
        f.write(outBedNew + '\n')


# Function for making a matrix instead of a bed file
def makematrix(region):
    global count, percent, totalLen, end_time, start_time, readName, forward, reverse, outBedNew, outBed, bufferPos, \
        bufferNeg, f, bedLen, isTabix

    # Split the input bed
    if args.debug:
        print(region)
    chrom = region.split('\t')[0]
    start = int(region.split('\t')[1])
    end = int(region.split('\t')[2])
    if args.debug:
        print(count)
    count += 1
    cnt = c.defaultdict(lambda: 0)
    cntneg = c.defaultdict(lambda: 0)
    totalreads = 0
    if args.smooth > 1:
        smooth = (int(args.smooth) / 2)
    else:
        smooth = 1

    # If a size is specified change the start and end sites to be Half way between +/- size, if not
    oldstart = start
    oldend = end
    if args.size > 10:
        oldstart = start
        oldend = end
        start = int((oldstart + (math.floor((oldend - oldstart) / 2))) - args.size)
        end = int((oldstart + (math.floor((oldend - oldstart) / 2))) + args.size)
    else:
        start = int((oldstart + (math.floor((oldend - oldstart) / 2))) - 10)
        end = int((oldstart + (math.floor((oldend - oldstart) / 2))) + 10)
    matrixline = []

    # On the first cycle, add a header line - Keep outside of main loop, otherwise it can fail on files with a bad first
    # region
    if count == 1:
        matrixline.append("Region")
        if args.centipede:
            for x in range(1, int(((end - start) + 1) * 2) + 1):
                matrixline.append(str(x))
        else:
            for x in range(0, int((end - start) + 1)):
                matrixline.append(str(x))
        matrixline.append("Total")
        f.write('\t'.join(map(str, matrixline)) + '\n')

    # Only process reads which are within the region boundaries
    if (start > 0 and smooth == 1) or (smooth > 1 and (start - smooth) > 0):

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
        # Pre-fill the counter with all base pairs and if centipede support is needed, make two counters
        # Centipede requires two matrices side by side, one with counts from positive strand reads and the other with
        # Negative strand reads.
        # if args.centipede:
        #     for i in range(int(start) - smooth, int(end + 1) + smooth):
        #         cnt[i] = 0
        #     for i in range(int(start) - smooth, int(end + 1) + smooth):
        #         cntneg[i] = 0
        # else:
        #     for i in range(int(start) - smooth, int(end + 1) + smooth):
        #         cnt[i] = 0
        matrixline = []

        # Try and get the region name, if there is no name, give it a number
        try:
            name = region.split('\t')[3]
        except IndexError:
            name = "region-" + str(count)

        # Try and extract strand information, if non existent or not '+'/'-', assume positive strand
        try:
            strand = region.split('\t')[5].rstrip('\n')
            if strand != '+' or strand != '-':
                strand = '+'
        except IndexError:
            strand = '+'

        # If we aren't doing a pileup, get the cut sites instead
        if not args.pileup:
            if isTabix:
                try:
                    for reg in datafile.fetch(chrom, start, end):
                        sreg = reg.split('\t')
                        try:
                            if sreg[5] == '+':
                                if sreg[0] == chrom and start <= int(sreg[1]) <= end:
                                    if args.binbed:
                                        cnt[int(sreg[1])] += 1
                                    try:
                                        cnt[int(sreg[1])] += float(sreg[4])
                                    except (IndexError, ValueError):
                                        cnt[int(sreg[1])] += 1
                            elif sreg[5] == '-':
                                if sreg[0] == chrom and start <= int(sreg[2]) <= end:
                                    if args.centipede:
                                        if args.binbed:
                                            cntneg[int(sreg[2])] += 1
                                        try:
                                            cntneg[int(sreg[2])] += float(sreg[4])
                                        except (IndexError, ValueError):
                                            cntneg[int(sreg[2])] += 1
                                    else:
                                        if args.binbed:
                                            cnt[int(sreg[2])] += 1
                                        try:
                                            cnt[int(sreg[2])] += float(sreg[4])
                                        except (IndexError, ValueError):
                                            cnt[int(sreg[2])] += 1
                            else:
                                raise IndexError
                        except IndexError:
                            if sreg[0] == chrom and start <= int(sreg[1]) <= end:
                                if args.binbed:
                                    cnt[int(sreg[1])] += 1
                                try:
                                    cnt[int(sreg[1])] += float(sreg[4])
                                except (IndexError, ValueError):
                                    cnt[int(sreg[1])] += 1
                except ValueError:
                    if args.debug:
                        print("Region" + chrom + " " + str(start) + " " + str(end) + " does not exist in this file")
            if not isTabix:
                for i in datafile.fetch(chrom, start + bufferPos - smooth, end - bufferNeg + smooth):
                    if not i.is_reverse:
                        forward += 1
                        if start <= int(i.pos + bufferPos) <= end:
                            cnt[i.pos + bufferPos] += 1
                    if i.is_reverse and args.centipede:
                        reverse += 1
                        if start <= int(i.pos + i.inferred_length + bufferNeg) <= end:
                            cntneg[i.pos + i.inferred_length + bufferNeg] += 1
                    if i.is_reverse and not args.centipede:
                        reverse += 1
                        if start <= int(i.pos + i.inferred_length + bufferNeg) <= end:
                            cnt[i.pos + i.inferred_length + bufferNeg] += 1

        # If we are doing a pileup, pileup!
        if args.pileup:
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

            # Support for SAM and BAM files
            if not isTabix:
                if strand == '+':
                    matrixline.append(name.rstrip('\n'))
                cols = 0
                nextcol = start
                if start > 0:
                    for pileupColumn in datafile.pileup(chrom, start, end):
                        if args.debug:
                            print(start, pileupColumn.pos, end)
                        if int(start) > int(pileupColumn.pos) or int(end) < int(pileupColumn.pos):
                            continue
                        while nextcol < pileupColumn.pos:
                            matrixline.append("0")
                            nextcol += 1
                            cols += 1
                        if args.debug:
                            print(pileupColumn.pos)
                        matrixline.append(pileupColumn.n)
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
                        matrixline.append("0")
                if strand == '-':
                    matrixline.append(name)
                    matrixline = matrixline[::-1]

        # If not doing pileup, append the counts, in order, to the line to be printed
        if not args.pileup or isTabix:
            # Support for centipede
            if args.centipede:
                for base in sorted(cnt.keys()):
                    if base < start or base > end:
                        continue
                    totalreads += int(cnt[base])
                    if smooth > 1:
                        smoothtot = 0
                        for i in range((base - smooth), (base + smooth)):
                            smoothtot += cnt[i]
                        matrixline.append(str(smoothtot / ((smooth * 2) + 1)))
                        continue
                    else:
                        matrixline.append(str(cnt[base]).rstrip('\n'))
                        continue
                revkeys = sorted(cntneg.keys())[::-1]
                for base in revkeys:
                    if base < start or base > end:
                        continue
                    totalreads += int(cntneg[base])
                    if smooth > 1:
                        smoothtot = 0
                        for i in range((base - smooth), (base + smooth)):
                            smoothtot += cntneg[i]
                        matrixline.append(str(smoothtot / ((smooth * 2) + 1)))
                    else:
                        matrixline.append(str(cntneg[base]).rstrip('\n'))

            if strand == '+' and not args.centipede:
                matrixline.append(name.rstrip('\n'))
                #for base in sorted(cnt.keys()):
                for base in range(start, end + 1, 1):
                    if base < start or base > end:
                        continue
                    totalreads += int(cnt[base])
                    if smooth > 1:
                        smoothtot = 0
                        for i in range((base - smooth), (base + smooth)):
                            smoothtot += cnt[i]
                        matrixline.append(str(smoothtot / ((smooth * 2) + 1)))
                    else:
                        matrixline.append(str(cnt[base]).rstrip('\n'))

            if strand == '-' and not args.centipede:
                #revkeys = sorted(cnt.keys())
                #for base in revkeys:
                for base in range(end + 1, start, 1):

                    if base < start or base > end:
                        continue
                    totalreads += int(cnt[base])
                    if smooth > 1:
                        smoothtot = 0
                        for i in range((base - smooth), (base + smooth)):
                            smoothtot += cnt[i]
                        matrixline.append(str(smoothtot / ((smooth * 2) + 1)))
                    else:
                        matrixline.append(str(cnt[base]).rstrip('\n'))
                matrixline.append(name)
                matrixline = matrixline[::-1]

        # Finish the line with a total count and print

        if args.sumsOnly:
            matrixline.append(totalreads)
            try:
                for i in range(1, len(matrixline)):
                    mSums[i] = mSums[i] + float(matrixline[i])
            except IndexError:
                for i in range(1, len(matrixline)):
                    mSums.append(float(matrixline[i]))
        else:
            matrixline.append(totalreads)
            if args.centipede:
                if len(matrixline) == ((args.size * 4) + 4):
                    f.write('\t'.join(map(str, matrixline)) + '\n')
                elif not args.quiet:
                    # TODO: Track down this bug!
                    print(str(region) + "Line was too short:- " + str(((args.size * 4) + 4))) + " > " + str(len(matrixline))
                    print("Please e-mail the command used and these files to kristofer.davie@kuleuven.be")
            else:
                if len(matrixline) == ((args.size * 2) + 3):
                    f.write('\t'.join(map(str, matrixline)) + '\n')
                elif not args.quiet:
                    print(str(region) + "Line was too short:- " + str(((args.size * 2) + 3))) + " > " + str(len(matrixline))
                    print("Please e-mail the command used and these files to kristofer.davie@kuleuven.be")
                    if args.debug:
                        print(matrixline)

if not args.quiet:
    print("Beginning to process file")

# Feedback and run the main loops
if not args.matrix:
    if not args.quiet:
        print("Counting Reads...")
    totalLen = (int(pysam.view("-c", args.inputFile)[0]))
    if not args.quiet:
        print("Total Number of Reads: " + str(totalLen))
    # Prevent error in PyCharm about unresolvable variable. Will exit if cannot open file.
    # noinspection PyUnboundLocalVariable
    for read in datafile:
        getsite(read)
    end_time = time.time()
    if not args.quiet:
        print("Processed 100% (" + str(count) + " of " + str(totalLen) + " reads) in " + str(
              end_time - start_time) + " seconds.")
elif args.matrix:
    for line in b:
        bedLen += 1
    b.seek(0)
    # Another PyCharm fix
    # noinspection PyUnboundLocalVariable
    for r in b:
        makematrix(r)

if args.sumsOnly:
    f.write('\t'.join(map(str, mSums)) + '\n')

    import matplotlib.pyplot as plt
    xT = []
    for x in range(-(int(args.size))-1, int(args.size), 1):
        xT.append(x)
    p = plt.plot(xT, mSums[1:-1], c='black')
    plt.ylabel('Aggregate reads')
    plt.xlabel('Position')
    plt.xlim((-(int(args.size)), (int(args.size))))
    plt.axvline(x=0, ls='--', lw=1)
    plt.savefig('.'.join([args.outputFile, 'pdf']), bbox_inches='tight')
f.close()
if not args.quiet and not args.matrix:
    print("Processing Finished")
    print("Total Forward Reads: " + str(forward))
    print("Total Reverse Reads: " + str(reverse))
    print("Reads not processed (unmapped etc): " + str(count - forward - reverse))
