#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from Bio import SeqIO
import numpy as np


def nucleotidefeature(fa, pos):
    '''
    This function is to calculate the single and double nucleotide content 
    for gene body region of the primary transcript of each annotated gene. 
    fa is the cds sequence of input primary transcript file.  
    '''
    fasta_sequences = SeqIO.parse(open(fa), 'fasta')
    mdict = {}
    n = 0
    xlist = ['Gene']
    two = []
    for i in ['A', 'T', 'C', 'G']:
        xlist.append(i + '-' + pos)

    for i in ['A', 'T', 'C', 'G']:
        for j in ['A', 'T', 'C', 'G']:
            string = i + j
            two.append(string)
            xlist.append(string + '-' + pos)
            n += 1

    for fasta in fasta_sequences:
        fname, sequence = fasta.id, str(fasta.seq)
        sequence = sequence.upper()
        if '.' in fname:
            n = fname.split('.')
            fname = n[0] + '.' + n[1]
        else:
            fname = fname
        name = fname.split('-')[0]
        tlist = []
        a = sequence.count('A')
        t = sequence.count('T')
        c = sequence.count('C')
        g = sequence.count('G')
        n = sequence.count('N')
        tlist = [a, t, c, g]
        for i in two:
            tlist.append(sequence.count(i))
        # ignoring uncertain nucleotide in the sequence
        tlist.append(len(sequence) - n)
        tlist = np.array(tlist)
        mdict[name] = tlist

    ndict = {}
    for i in sorted(mdict):
        if i not in ndict:
            ndict[i] = []
        for j in mdict[i][:-1]:
            ndict[i].append(str(j / float(mdict[i][-1])))

    fasta_sequences.close()
    # return the title and dictionary with gene name and appended values.
    return xlist, ndict

'''
def proteinfeature(protein):

    This function is to calculate the amino acid content for each gene, protein
    is the amino acid sequence of primary transcript
    
    fasta_sequences = SeqIO.parse(open(protein), 'fasta')

    mdict = {}
    num = 0
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        if '.' in name:
            n = name.split('.')
            name = n[0] + '.' + n[1]
        else:
            name = name

        a = sequence.count('A')
        r = sequence.count('R')
        n = sequence.count('N')
        d = sequence.count('D')
        c = sequence.count('C')
        q = sequence.count('Q')
        e = sequence.count('E')
        g = sequence.count('G')
        h = sequence.count('H')
        i = sequence.count('I')
        l = sequence.count('L')
        k = sequence.count('K')
        m = sequence.count('M')
        f = sequence.count('F')
        p = sequence.count('P')
        s = sequence.count('S')
        t = sequence.count('T')
        w = sequence.count('W')
        y = sequence.count('Y')
        v = sequence.count('V')
        mdict[name] = np.array([a, r, n, d, c, q, e, g, h, i, l, k, m, f, p, s, t, w, y, v, len(sequence)])

    ndict = {}
    xlist = ['Gene', 'Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile',
             'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val',
             'Alaratio', 'Argratio', 'Asnratio', 'Aspratio', 'Cysratio', 'Glnratio', 'Gluratio',
             'Glyratio', 'Hisratio', 'Ileratio', 'Leuratio', 'Lysratio', 'Metratio', 'Pheratio',
             'Proratio', 'Serratio', 'Thrratio', 'Trpratio', 'Tyrratio', 'Valratio']
    for i in sorted(mdict):
        if i not in ndict:
            ndict[i] = []
        for j in mdict[i][:-1]:
            ndict[i].append(str(j))
        for j in mdict[i][:-1]:
            ndict[i].append(str(j / float(mdict[i][-1])))
    num += 1

    fasta_sequences.close()
    return xlist, ndict
'''

def merge(*args):
    '''
    Merge multiple features of each gene together. Each argument stands for the
    dictionary of the class of each feature for each gene
    '''
    mdict = {}

    # fine the common keys in multiple dictionaries
    common_keys = set(args[0].keys())
    for d in args[1:]:
        common_keys.intersection_update(set(d.keys()))

    # appending values through common_keys to a dictionary
    for i in common_keys:
        if i not in mdict:
            mdict[i] = []
        for d in args[:]:
            mdict[i].extend(d.get(i))

    return mdict
