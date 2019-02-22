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
            with open('seqFile1','w+') as out:
                out.write(seq1)
            with open('seqFile2','w+') as out2:
                out2.write(seq2)
            realignCommand = "%s | cPecanRealign %s %s %s" % (cig, "--rescoreByPosteriorProb", 'seqFile1', 'seqFile2')
            for realignLineByPosteriorProb in [i for i in popenCatch(realignCommand).split("\n") if i != '']:
                line = cigarReadFromString(realignLineByPosteriorProb)
                print(line)
    return cigars


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
                print(int(pos[2]) + int(pos[3]))
                print(start)
                print(int(pos[2]))
                if int(pos[2]) + int(pos[3]) >= start >= int(pos[2]):
                    allBlocks.append(current)
                elif int(pos[2]) >= start and int(pos[3]) <= end:
                    allBlocks.append(current)
                elif int(pos[2]) <= end <= int(pos[2]) + int(pos[3]) or int(pos[2]) <= start and int(pos[3]) >= end:
                    allBlocks.append(current)
                    break
            elif inBlock:
                current.append(line)
    for i in range(len(allBlocks)):
        temp = []
        for line in allBlocks[i]:
            if line[0] != '#' and line[0] != 'a':
                temp.append(line.split())
        allBlocks[i] = temp
    if len(allBlocks) != 1:
        first = allBlocks[0]
        last = allBlocks[len(allBlocks)-1]
        for i in range(len(first)):
            allBlocks[0][i][6] = first[i][6][start - int(first[i][2]):]
        for i in range(len(last)):
            allBlocks[len(allBlocks)-1][i][6] = last[i][6][:int(last[i][3])+int(last[i][2])-end]
    else:
        block = allBlocks[0]
        for i in range(len(block)):
            allBlocks[0][i][6] = block[i][6][start-int(block[2]):int(block[3])+int(block[2])-end]
    return allBlocks


def cigar(seq1, seq2):
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
            full += 'I'
        elif seq1[i] == '-' and '-' != seq2[i]:
            full += 'D'
    count = 1
    cigar = ''
    for i in range(1,len(seq1)):
        if full[i-1] == full[i]:
            count += 1
        else:
            cigar += str(count) + full[i-1]
            count = 1
    cigar += str(count) + full[len(seq1)]
    return cigar


def parseArgs():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-p',help='region of interest, e.g. 400-500')
    parser.add_argument('-m',help='maf file')
    parser.add_argument('-r',help='reference genome as written in maf file')
    parser.add_argument('-t',help='target genomes as written in maf file')
    #parser.add_argument('-genomes',help='genome directory')
    #parser.add_argument('--blockSize',type=int,default=10000,help='size of upstream and downstream for HMM to train')
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
    print(dictBlock)
    pairwise(opts.r,dictBlock,genomes)



if __name__ == '__main__':
    main()
