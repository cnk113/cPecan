#!/usr/bin/env python3
# Chang Kim

'''
Pulls overlapping MAF blocks from alignments
Trains the pair HMM on the alignment block and surrounding sequence
'''

import argparse
import os
#from sonLib.bioio import popenCatch, cigarReadFromString


def pairwise(ref,blocks,targets):
    '''
    :return:
    '''
    cigars = []
    for block in blocks:
        seq1 = block.get(ref)
        for target in targets:
            seq2 = block.get(target)
            cig = cigar(seq1, seq2)
            realignCommand = "%s | cPecanRealign %s %s %s" % (cig, "--outputPosteriorProbs", 'seqFile1', 'seqFile2')
            for realignLineByPosteriorProb in [i for i in popenCatch(realignCommand).split("\n") if i != '']:
                line = cigarReadFromString(realignLineByPosteriorProb)
                print(line)
    return cigars


def cigar(seq1, seq2, ran):
    '''
    :param seq1: target
    :param seq2: query
    :return: cigar
    '''
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    full = ''
    for i in range(len(seq1)):
        if seq1[i] != '-' and '-' != seq2[i]:
            full += 'M'
        elif seq1[i] != '-' and seq2[i] == '-':
            full += 'D'
        elif seq1[i] == '-' and '-' != seq2[i]:
            full += 'I'
    count = 1
    r = ran.split('-')
    cigar = 'cigar: query' + r[0] + ' ' + r[1]
    for i in range(1,len(seq1)):
        if full[i-1] == full[i]:
            count += 1
        else:
            cigar += full[i-1] + ' ' + str(count) + ' '
            count = 1
    cigar += full[len(seq1)] + ' ' + str(count)
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
            allBlocks[0][i][6] = first[i][6][startOffset:]
        for i in range(len(last)):
            if endOffset == 0:
                allBlocks[len(allBlocks) - 1][i][6] = last[i][6]
            else:
                allBlocks[len(allBlocks)-1][i][6] = last[i][6][:endOffset]
    else:
        block = allBlocks[0]
        for i in range(len(block)):
            if endOffset == 0:
                allBlocks[0][i][6] = block[i][6][startOffset:]
            else:
                allBlocks[0][i][6] = block[i][6][startOffset:endOffset]
    return allBlocks


def parseArgs():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-p',help='region of interest, e.g. 400-500')
    parser.add_argument('-m',help='maf file')
    parser.add_argument('-r',help='reference genome as written in maf file')
    parser.add_argument('-t',help='target genomes as written in maf file')
    return parser.parse_args()


def main():
    opts = parseArgs()
    with open(opts.t,'r') as infile:
        genomes = [line.rstrip() for line in infile]
    blocks = parseMaf(opts.m, opts.p)
    dictBlock = []
    for block in blocks:
        map = {}
        for line in block:
            map[line[1].split('.')[0]] = line[6]
        dictBlock.append(map)
    #print(dictBlock)
    pairwise(opts.r,dictBlock,genomes)



if __name__ == '__main__':
    main()
