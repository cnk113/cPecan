#!/usr/bin/env python3
# Chang Kim

'''
Pulls overlapping MAF blocks from alignments
Trains the pair HMM on the alignment block and surrounding sequence
'''

import argparse
import os
from sonLib.bioio import popenCatch

class MafBlock:
    '''
    Contains the MAF block in a dictionary where the species is the key
    position function is relative to the reference genome position
    '''
    def __init__(self,range,):
        '''
        just instantiates the block
        :param range:
        '''
        self.genomes = {}
        self.range = range

    def createBlock(self,maf):
        '''
        actually instantiates the block, since we don't want to do it for all MAF blocks in the MAF file
        :param maf:
        :return:
        '''
        lines = maf[self.range[0]:self.range[1]]
        for line in lines:
            attr = line.split()
            genome = attr[1].split('.')
            self.genomes[genome[0]] = (genome[1],attr[2],attr[3],attr[4],attr[6].rstrip())

    def getReference(self,genomes,size):
        '''
        gets the surrounding sequences of the alignment blocks
        :param genomes:
        :param size:
        :return:
        '''
        for genomes in genomes:



    def get(self):
        '''
        :return: the alignment block
        '''
        return self.genomes


def parseMaf(maf,range):
    '''
    :param maf: maf file
    :param range: search range
    :return: list of MafBlocks in sorted order of the reference genome position
    '''
    r = range.split('-')
    start = int(r[0])
    end = int(r[1])
    with open(maf) as infile:
        inBlock = False
        allBlocks = []
        current = []
        for line in infile:
            if line[0] == 'a' and not inBlock:
                current.append(line)
                inBlock = True
            elif not line and inBlock:
                inBlock = False
                pos = current[1].split()
                if int(pos[2]) + int(pos[3]) >= start >= int(pos[2]):
                    allBlocks.append(current)
                elif int(pos[2]) >= start and int(pos[3]) <= end:
                    allBlocks.append(current)
                elif int(pos[2]) <= end <= int(pos[2]) + int(pos[3]):
                    allBlocks.append(current)
                    break
            elif inBlock:
                current.append(line)
    mafBlocks = []
    for block in allBlocks:
        mafBlocks.append(MafBlock(block[1:]))
    return mafBlocks


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
        elif seq1[i] == '-' and '-' == seq2[i]:
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
    parser.add_argument('-pos',help='region of interest, e.g. 400-500')
    parser.add_argument('-maf',help='maf file')
    parser.add_argument('-genomes',help='genome directory')
    #parser.add_argument('--blockSize',type=int,default=10000,help='size of upstream and downstream for HMM to train')
    return parser.parse_args()


def main():
    opts = parseArgs()
    blocks = parseMaf(opts.maf,opts.pos)
    genomes = {}
    for f in os.listdir(opts.genomes):
        if os.path.isfile(os.path.join(opts.genomes,f)):
            genomes[f.split('.')[0]] = Fasta(os.path.join(opts.genomes,f))
    for block in blocks:
        block.getReference(genomes,opts.blockSize)

    cig = cigar(seq1,seq2)
    realignCommand = "%s | cPecanRealign %s %s %s" % (cig, "--rescoreByPosteriorProb", seqFile1, seqFile2)
    output = [i for i in popenCatch(realignCommand).split("\n") if i != '']
    realignCigar = cigarReadFromString(output)



if __name__ == '__main__':
    main()
