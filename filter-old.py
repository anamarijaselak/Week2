import sys
import os
import argparse
import math
import time
import os.path

import pylab as pl
import numpy as np
from Bio import SeqIO
from _collections import defaultdict


def check_arguments():
    parser = argparse.ArgumentParser(usage='Script for analysis of fastq files')
    if(len(sys.argv) != 3):
        parser.error('''Invalid number of arguments!
            Needed arguments : <INPUT-FA-FILE> <OUTPUT-FOLDER>
            INPUT-FA-FILE - fasta file that is going to be processed
            OUTPUT-FOLDER - directory in which the results are going to be stored''')
        exit()
    parser.add_argument("fasta", help="Fasta file that needs to be processed")
    parser.add_argument("outputDir")
    args = parser.parse_args()
    if not (args.fasta.endswith('.fasta') or args.fasta.endswith('.fa')) or (not os.path.isdir(args.outputDir)):
        parser.error('''Invalid types of given arguments!
            Needed arguments : <INPUT-FA-FILE> <OUTPUT-FOLDER>
            INPUT-FA-FILE - fasta file that is going to be processed
            OUTPUT-FOLDER - directory in which the results are going to be stored''')
        exit()


def makeHistogram(destinationFolder, entropiesDict):
    pl.hist(entropiesDict.keys(), len(entropiesDict), weights=entropiesDict.values())
    pl.title("Sequence entropy distribution")
    pl.savefig(destinationFolder + '/FastaEntropiesHist.png')


def human_readable(bytes):
    size = ""
    size = str(bytes / 1000000000.)
    return size + "\n"


def statistics(fastaPath, destinationFile, numberOfReads, numberOfBadReads, numberOfKeptReads, start_time):
    fastaFile = open(fastaPath)
    statisticsFile = open(destinationFile + "/FastaTimeStatistics.txt", "a+")
    bytes = os.path.getsize(fastaPath)
    size = human_readable(bytes)
    statisticsFile.write(str(size))
    statisticsFile.write(str(time.time()-start_time)+"\n")


def main():
    total_entropy_calc_time = 0.
    total_dict_time = 0.
    start_time=time.time()
    check_arguments()
    entropiesDict = defaultdict(int)
    max_entropy = 0.
    numberOfReads = 0
    numberOfKeptReads = 0
    numberOfBadReads = 0
    newFastaFile = open(sys.argv[2] + "/newFastafile.fa", "w+")
    for record in SeqIO.parse(sys.argv[1], "fasta"):
        numberOfReads = numberOfReads + 1
        basesInRecord = 0
        baseDict = defaultdict(int)
        if "N" in record.seq :
            numberOfBadReads = numberOfBadReads + 1
            continue
        entropy = 0
        for base in record.seq:
            baseDict[base] = baseDict[base] + 1
            basesInRecord = basesInRecord + 1
        for key in baseDict.keys():
            start = time.time()
            entropy = entropy + (baseDict[key] / float(basesInRecord)) * math.log(baseDict[key] / float(basesInRecord), 2)
            stop = time.time()
            total_entropy_calc_time += (stop-start)
            start = time.time()
            entropiesDict[round(-1 * entropy, 2)] = entropiesDict[round(-1 * entropy, 2)] + 1
            stop = time.time()
            total_dict_time += (stop-start)
        if -1 * entropy > 0.5:
            newFastaFile.write(">")
            newFastaFile.write(record.name + "\n")
            newFastaFile.write(str(record.seq) + "\n")
            numberOfKeptReads = numberOfKeptReads + 1
        else :
            numberOfBadReads = numberOfBadReads + 1
    makeHistogram(sys.argv[2], entropiesDict)
    statistics(sys.argv[1], sys.argv[2], numberOfReads, numberOfBadReads, numberOfKeptReads, start_time)
    print 'Entropy calculation', total_entropy_calc_time
    print 'Dictionary access', total_dict_time

if __name__ == '__main__':
    main()
