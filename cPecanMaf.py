#!/usr/bin/env python3
# Chang Kim

'''
Pulls overlapping MAF blocks from alignments
Trains the pair HMM on the alignment block and surrounding sequence
'''

import argparse
import os
from sonLib.bioio import popenCatch, cigarReadFromString


def pairwise(ref,blocks,targets,ran):
    '''
    :return:
    '''
    cigars = []
    for block in blocks:
        seq1 = block.get(ref)
        seqFile1 = targets.get(ref)
        for species in block:
            if species != ref:
                seq2 = block.get(species)
                seqFile2 = targets.get(species)
                cig = cigar(seq1,seq2,ran)
                realignCommand = "%s | cPecanRealign %s %s %s" % (cig, "--outputPosteriorProbs", seqFile1, seqFile2)
                for realignLineByPosteriorProb in [i for i in popenCatch(realignCommand).split("\n") if i != '']:
                    line = cigarReadFromString(realignLineByPosteriorProb)
                    print(line)
    return cigars


def cigar(seq1, seq2):
    '''
    :param seq1: ref
    :param seq2: query
    :return: cigar
    '''
    cigar = 'cigar: query ' + seq2[2] + ' ' + str(seq2[0]) + ' ' + str(seq2[0] + seq2[1] + 1) + ' '
    seq1 = seq1[3].upper()
    seq2 = seq2[3].upper()
    full = ''
    for i in range(len(seq1)):
        if seq1[i] != '-' and '-' != seq2[i]:
            full += 'M'
        elif seq1[i] != '-' and seq2[i] == '-':
            full += 'D'
        elif seq1[i] == '-' and '-' != seq2[i]:
            full += 'I'
    count = 1
    for i in range(1,len(seq1)):
        if full[i-1] == full[i]:
            count += 1
        else:
            cigar += full[i-1] + ' ' + str(count) + ' '
            count = 1
    cigar += full[len(seq1)-1] + ' ' + str(count)
    return cigar


def parseMaf(maf,ran):
    '''
    :param maf: maf file
    :param ran: search range
    :return: list of MafBlocks in sorted order of the reference genome position
    '''
    r = ran.split('-')
    start = int(r[0])
    end = int(r[1])
    with open(maf) as infile:
        inBlock = False
        allBlocks = []
        current = []
        for line in infile:
            line = line.rstrip()
            if line and line[0] == 'a' and not inBlock:
                current.append(line)
                inBlock = True
            elif not line and inBlock:
                inBlock = False
                pos = current[1].split()
                a = int(pos[2])
                b = a + int(pos[3])
                if a <= start < b < end:
                    allBlocks.append(current)
                elif start <= a < b < end:
                    allBlocks.append(current)
                elif start <= a < end <= b or a <= start < end <= b:
                    allBlocks.append(current)
                    break
                current = []
            elif inBlock:
                current.append(line)
    for i in range(len(allBlocks)):
        temp = []
        for line in allBlocks[i]:
            if line[0] != '#' and line[0] != 'a':
                temp.append(line.split())
        allBlocks[i] = temp

    startOffset = start - int(allBlocks[0][0][2])
    endOffset = end - (int(allBlocks[len(allBlocks) - 1][0][3]) + int(allBlocks[len(allBlocks) - 1][0][2]))
    if len(allBlocks) != 1:
        first = allBlocks[0]
        last = allBlocks[len(allBlocks)-1]
        for i in range(len(first)):
            allBlocks[0][i][6] = first[i][6][startOffset:] # Cutting off alignment to our query range
            allBlocks[0][i][2] = int(allBlocks[0][i][2]) + startOffset # adjust start position of alignment in fasta
        for i in range(len(last)):
            if endOffset == 0:
                allBlocks[len(allBlocks) - 1][i][6] = last[i][6]
            else:
                allBlocks[len(allBlocks)-1][i][6] = last[i][6][:endOffset]
            allBlocks[0][i][3] = int(allBlocks[0][i][3]) - abs(endOffset) # adjust length of alignment
    else:
        block = allBlocks[0]
        for i in range(len(block)):
            if endOffset == 0:
                allBlocks[0][i][6] = block[i][6][startOffset:]
            else:
                allBlocks[0][i][6] = block[i][6][startOffset:endOffset]
            allBlocks[0][i][2] = int(allBlocks[0][i][2]) + startOffset
            allBlocks[0][i][3] = int(allBlocks[0][i][3]) - abs(endOffset)
    return allBlocks


def parseArgs():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-p',help='region of interest, e.g. 400-500')
    parser.add_argument('-m',help='maf file')
    parser.add_argument('-r',help='reference genome as written in maf file')
    parser.add_argument('-t',help='target genomes as written in maf file and their corresponding references')
    return parser.parse_args()


def main():
    opts = parseArgs()
    with open(opts.t,'r') as infile:
        ref = {}
        for line in infile:
            attr = line.rstrip().split()
            ref[attr[0]] = attr[1]
    blocks = parseMaf(opts.m, opts.p)
    dictBlock = []
    for block in blocks:
        map = {}
        for line in block:
            print line
            map[line[1].split('.')[0]] = (line[2],line[3],line[4],line[6])
        dictBlock.append(map)
    #print(dictBlock)
    pairwise(opts.r,dictBlock,ref,opts.p)



if __name__ == '__main__':
    main()
