#!/usr/bin/env python3
# Chang Kim

'''
Pulls overlapping MAF blocks from alignments
Trains the pair HMM on the alignment block and surrounding sequence
'''

import argparse
import os
from sonLib.bioio import popenCatch, cigarReadFromString, system


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
                if seq2[2] == '-':
                    continue


                seqFile2 = targets.get(species)
                cig = cigar(seq1,seq2)
                #realignCommand = "echo '%s' | /mnt/c/Users/Chang/Desktop/allcode/sonLib/bin/cPecanRealign %s %s %s" % (cig, "-u /dev/stdout", seqFile1, seqFile2)

                print popenCatch(realignCommand).split("\n")
                #for realignLineByPosteriorProb in [i for i in popenCatch(realignCommand).split("\n") if i != '']:
                    #line = cigarReadFromString(realignLineByPosteriorProb)
                    #print(line)
    return cigars

'''
def cigar(seq1, seq2):
    print seq1
    print seq2
    cigar = 'cigar: ' + seq1[4] +' '+ str(seq1[0]) + ' ' + str(seq1[1]) + ' ' + seq1[2] + ' ' +seq2[4]+ ' ' + str(seq2[0]) + ' ' + str(seq2[1]) + ' ' + seq2[2] + ' 0 '
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
'''
def psl2cigar(files):
    '''
    '''
    cigars = {}
    for psl in files:
        with open(psl) as infile:
            for line in psl:
                attr = line.rstrip().split()
                
    return cigars


def axtChain(axt, fasta, ref):
    '''
    :param axt:
    :param fasta:
    :param ref:
    '''
    for aln in axt:
        cmd = 'axtChain -linearGap=medium ' + axt.get(aln) + ' -faQ ' + fasta.get(ref) + ' -faT ' + fasta.get(aln) + aln + '.chain'
        print cmd
        system(cmd)
    for aln in axt:
        cmd = 'chain2psl ' + aln + '.chain'
        

def maf2axt(maf, genomes):
    '''
    :param maf: maf file
    :param genomes: genomes to align
    '''
    files = {}
    for genome in genomes:
        files[genome] = genome + '.axt'
    with open(maf) as infile:
        newBlock = True
        for line in infile:
            line = line.rstrip()
            if line:
                if newBlock and line[0] == 's':
                    currentRef = line.split()
                    refSpecies = currentRef[1].split('.')
                    newBlock = False
                elif line[0] == 's' and not newBlock:
                    target = line.split()
                    species = target[1].split('.')
                    if species[0][0:3] != 'Anc':
                        with open(files.get(species[0]), 'a') as out:
                            header = '0 ' + '.'.join(refSpecies[1:]) + ' ' + currentRef[2] + ' ' + str(int(currentRef[2]) + int(currentRef[3])) + ' '
                            header += '.'.join(species[1:]) + ' ' + target[2] + ' ' + str(int(target[2]) + int(target[3])) + ' ' + target[4] + ' ' + '1000\n'
                            out.write(header + currentRef[6] + '\n' + target[6] + '\n\n')
            else:
                newBlock = True
    return files, refSpecies[0]


def parseArgs():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-m',help='maf file')
    parser.add_argument('-t',help='target genomes as written in maf file and their corresponding references')
    return parser.parse_args()


def main():
    opts = parseArgs()
    with open(opts.t,'r') as infile:
        ref = {}
        for line in infile:
            attr = line.rstrip().split()
            ref[attr[0]] = attr[1]
    axt, refSpecies = maf2axt(opts.m, ref)
    axtChain(axt, ref, refSpecies)
    psl2cigar(axt)


if __name__ == '__main__':
    main()
