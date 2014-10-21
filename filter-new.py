import os
import os.path
import sys
import math
import time
import argparse
from array import *
from collections import defaultdict

import pylab as pl
import numpy as np
from scipy import weave
from Bio import SeqIO


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
    if not (args.fasta.endswith('.fasta') or args.fasta.endswith('.fa') or args.fasta.endswith('fastq')) or (not os.path.isdir(args.outputDir)):
        parser.error('''Invalid types of given arguments!
            Needed arguments : <INPUT-FA-FILE> <OUTPUT-FOLDER>
            INPUT-FA-FILE - fasta file that is going to be processed
            OUTPUT-FOLDER - directory in which the results are going to be stored''')
        exit()


def makeHistogram(destinationFolder, entropiesDict):
    axes = pl.gca()
    bins = np.arange(200)
    width = 0.35
    pl.bar(bins, entropiesDict, width) # ovu bins varijablu sam mjenjala na 1000 i jedan nacin
    axes.set_xticks(bins+width)
    axes.set_xticklabels(['%.1f' % (n/1.) for n in bins])
    pl.title("Sequence entropy distribution")
    pl.savefig(destinationFolder + '/FastaEntropiesHist.png')


def human_readable(bytes):
    size = ""
    size = str(bytes / 1000000000.)
    return size + "\n"


def statistics(fastaPath, destinationFile, num_reads, numberOfBadReads, numberOfKeptReads, start_time):
    fastaFile = open(fastaPath)
    statisticsFile = open(destinationFile + "/FastaTimeStatistics.txt", "a+")
    bytes = os.path.getsize(fastaPath)
    size = human_readable(bytes)
    statisticsFile.write(str(size))
    statisticsFile.write(str(time.time()-start_time)+"\n")

def nucl_entropy(nucl_seq):
    basesInRecord = 0
    baseArray = array('f', [0,0,0,0])
    entropy = 0
    for base in nucl_seq:
        if(base=="A"):
            baseArray[0] += 1
        if(base=="C"):
            baseArray[1] += 1
        if(base=="G"):
            baseArray[2] += 1
        if(base=="T"):
            baseArray[3] += 1
        basesInRecord = basesInRecord + 1
    for i in range(len(baseArray)):
        prob = baseArray[i] / float(basesInRecord)
        entropy = entropy + (prob * math.log(prob, 2))
    return entropy


def main():
    total_entropy_calc_time = 0.
    total_dict_time = 0.
    start_time=time.time()
    check_arguments()
    entropies = array('f', [0]*200)
    max_entropy = 0.
    num_reads = 0
    numberOfKeptReads = 0
    numberOfBadReads = 0
    newFastaFile = open(sys.argv[2] + "/newFastafile.fa", "w+")
    for record in SeqIO.parse(sys.argv[1], "fasta"):
        if "N" in record.seq :
            numberOfBadReads = numberOfBadReads + 1
            continue
        entropy = nucl_entropy(str(record.seq))
        position = int(round(-1 * entropy, 2)*100)
        entropies[position] = entropies[position] + 1
        num_reads = num_reads + 1
        if -1 * entropy > 0.5:
            newFastaFile.write(">")
            newFastaFile.write(record.name + "\n")
            newFastaFile.write(str(record.seq) + "\n")
            numberOfKeptReads = numberOfKeptReads + 1
        else :
            numberOfBadReads = numberOfBadReads + 1

    makeHistogram(sys.argv[2], entropies)
    statistics(sys.argv[1], sys.argv[2], num_reads, numberOfBadReads, numberOfKeptReads, start_time)

if __name__ == '__main__':
    main()
