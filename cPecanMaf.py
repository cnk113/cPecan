# Chang Kim

import argparse
import os
from sonLib.bioio import popenCatch, cigarReadFromString, system


def realign(fasta, ref):
    '''
    :param fasta: genome fasta files
    :param ref: reference species
    '''
    seqFile1 = fasta.get(ref)
    for genome in fasta:
        if genome == ref:
            continue
        axt = genome + '.chained.axt'
        with open(axt) as infile:
            block = []
            for line in infile:
                if line:
                    col = line.rstrip().split()
                    block.append(col)
                    if len(block) == 3:
                        cig = axtToCigar(block)
                        seqFile2 = fasta.get(genome)
                        realignCommand = "echo '%s' | cPecanRealign %s %s %s" % (cig, "-u /dev/stdout", seqFile1, seqFile2)
                        print popenCatch(realignCommand).split("\n")
                else:
                    block = []


def axtToCigar(block):
    '''
    :param block: axt alignment block
    :return: cigar string
    '''
    header = block[0]
    cigar = 'cigar: ' + header[1] + ' ' + header[2] + ' ' + header[3] + ' + ' + header[4] + ' ' + header[5] + ' ' + header[6] + ' ' + header[7] + ' 0 '
    seq1 = block[1][0].upper()
    seq2 = block[2][0].upper()
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


def chain(axt, fasta, ref):
    '''
    :param axt: axt files
    :param fasta: genome reference files
    :param ref: reference species
    '''
    system('faToTwoBit ' + fasta.get(ref) + ' ' + aln + '.2bit')
    for aln in axt:
        system('axtChain -linearGap=medium ' + axt.get(aln) + ' -faQ ' + fasta.get(ref) + ' -faT ' + fasta.get(aln) + ' ' + aln + '.chain')
        system('faToTwoBit ' + fasta.get(aln) + ' ' + aln + '.2bit')
        system('chainToAxt ' + aln + '.chain ' + aln + '.2bit ' + ref + '.2bit ' + aln + '.chained.axt'


def mafToAxt(maf, genomes):
    '''
    :param maf: maf file
    :param genomes: genomes to align
    '''
    files = []
    idx = {}
    filedir = {}
    for genome in genomes:
        files.append(genome + '.axt')
        filedir[genome] = genome + '.axt'
        idx[genome] = 0
    with open(maf) as infile:
        filedata = {filename: open(filename,'w+') for filename in files}
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
                    if filedir.get(species[0]) != None:
                        header = str(idx.get(species[0])) + ' ' + '.'.join(refSpecies[1:]) + ' ' + currentRef[2] + ' ' + str(int(currentRef[2]) + int(currentRef[3])) + ' '
                        header += '.'.join(species[1:]) + ' ' + target[2] + ' ' + str(int(target[2]) + int(target[3])) + ' ' + target[4] + ' ' + '1000\n'
                        filedata.get(filedir.get(species[0])).write(header + currentRef[6] + '\n' + target[6] + '\n\n')
                        idx[species[0]] += 1
            else:
                newBlock = True
        for f in filedata.values():
            f.close()
    return filedir, refSpecies[0]


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
    axt, refSpecies = mafToAxt(opts.m, ref)
    chain(axt, ref, refSpecies)
    realign(ref, refSpecies)


if __name__ == '__main__':
    main()
