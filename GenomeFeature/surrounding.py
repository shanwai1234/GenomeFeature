#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from pyfasta import Fasta


def upstream(fa, ann, kb1, kb2):
    '''
    Extracting gene upstream sequences. fa is genome assembly file, ann is
    the annotation file, kb1 is the defined length of 5' UTR, kb2 is the defined
    length of upstream.
    '''
    f = Fasta(fa)
    fh = open(ann, 'r')
    out1 = open('5-UTR.fa', 'w')
    out2 = open('upstream.fa', 'w')
    mdict = {}
    ndict = {}
    for line in fh:
        # this is the demo line that we want to filter out
        # chr7    GLEAN   Gene    25420153        25421713        0.953889        -       .       Name=Pgl_GLEAN_10006696;
        if line.startswith('#'):
            continue
        new = line.strip().split('\t')
        if new[2] != 'CDS':
            continue
        n = new[-1].split(';')
        for i, j in enumerate(n):
            if 'Parent=' in j:
                mindex = i
        g = n[mindex].split('.')
        t = g[0].replace('Parent=', '')
        if '_' in t:
            gene = t.split('_')[0]
        else:
            gene = t
        if gene not in mdict:
            mdict[gene] = []
            ndict[gene] = []
        ndict[gene].append(new[0])
        ndict[gene].append(new[6])
        mdict[gene].append(int(new[3]))
        mdict[gene].append(int(new[4]))

    for gene in sorted(mdict):
        # when gene is on the positive strand
        if ndict[gene][1] == '+':    
            start = min(mdict[gene])
            start1 = start - (int(kb1) + 1)
            stop1 = start - 1
            start2 = start1 - (int(kb2) + 1)
            stop2 = start1 - 1
            k1 = f.sequence({'chr': ndict[gene][0], 'start': start1, 'stop': stop1, 'strand': ndict[gene][1]})
            out1.write('>{0}-5UTR'.format(gene) + '\n')
            out1.write(k1 + '\n')
            k2 = f.sequence({'chr': ndict[gene][0], 'start': start2, 'stop': stop2, 'strand': ndict[gene][1]})
            out2.write('>{0}-upstream'.format(gene) + '\n')
            out2.write(k2 + '\n')
        # when gene is on the negative strand
        elif ndict[gene][1] == '-':
            stop = max(mdict[gene])
            start1 = stop + 1
            stop1 = stop + (int(kb1) + 1)
            start2 = stop1 + 1
            stop2 = stop1 + (int(kb2) + 1)
            k1 = f.sequence({'chr': ndict[gene][0], 'start': start1, 'stop': stop1, 'strand': ndict[gene][1]})
            out1.write('>{0}-5UTR'.format(gene) + '\n')
            out1.write(k1 + '\n')
            k2 = f.sequence({'chr': ndict[gene][0], 'start': start2, 'stop': stop2, 'strand': ndict[gene][1]})
            out2.write('>{0}-upstream'.format(gene) + '\n')
            out2.write(k2 + '\n')
    fh.close()
    out1.close()
    out2.close()


def downstream(fa, ann, kb1, kb2):
    '''
    Extracting gene upstream sequences. fa is genome assembly file, ann is
    the annotation file, kb1 is the defined length of 3' UTR, kb2 is the defined
    length of downstream.
    '''
    f = Fasta(fa)
    fh = open(ann, 'r')
    out1 = open('3-UTR.fa', 'w')
    out2 = open('downstream.fa', 'w')
    mdict = {}
    ndict = {}
    for line in fh:
        # this is the demo line that we want to filter out
        # chr7    GLEAN   Gene    25420153        25421713        0.953889        -       .       Name=Pgl_GLEAN_10006696;
        if line.startswith('#'):
            continue
        new = line.strip().split('\t')
        if new[2] != 'CDS':
            continue
        n = new[-1].split(';')
        for i, j in enumerate(n):
            if 'Parent=' in j:
                mindex = i
        g = n[mindex].split('.')
        t = g[0].replace('Parent=', '')
        if '_' in t:
            gene = t.split('_')[0]
        else:
            gene = t
        if gene not in mdict:
            mdict[gene] = []
            ndict[gene] = []
        ndict[gene].append(new[0])
        ndict[gene].append(new[6])
        mdict[gene].append(int(new[3]))
        mdict[gene].append(int(new[4]))

    for gene in sorted(mdict):
        if ndict[gene][1] == '+':
            stop = max(mdict[gene])
            start1 = stop + 1
            stop1 = stop + (int(kb1) + 1)
            start2 = stop1 + 1
            stop2 = stop1 + (int(kb2) + 1)
            k1 = f.sequence({'chr': ndict[gene][0], 'start': start1, 'stop': stop1, 'strand': ndict[gene][1]})
            out1.write('>{0}-3UTR'.format(gene) + '\n')
            out1.write(k1 + '\n')
            k2 = f.sequence({'chr': ndict[gene][0], 'start': start2, 'stop': stop2, 'strand': ndict[gene][1]})
            out2.write('>{0}-downstream'.format(gene) + '\n')
            out2.write(k2 + '\n')
        elif ndict[gene][1] == '-':
            start = min(mdict[gene])
            start1 = start - (int(kb1) + 1)
            stop1 = start - 1
            start2 = start1 - (int(kb2) + 1)
            stop2 = start1 - 1
            k1 = f.sequence({'chr': ndict[gene][0], 'start': start1, 'stop': stop1, 'strand': ndict[gene][1]})
            out1.write('>{0}-3UTR'.format(gene) + '\n')
            out1.write(k1 + '\n')
            k2 = f.sequence({'chr': ndict[gene][0], 'start': start2, 'stop': stop2, 'strand': ndict[gene][1]})
            out2.write('>{0}-downstream'.format(gene) + '\n')
            out2.write(k2 + '\n')
    fh.close()
    out1.close()
    out2.close()


def cds(fa, ann):
    f = Fasta(fa)
    fh = open(ann, 'r')
    out1 = open('exon.fa', 'w')
    mdict = {}
    for line in fh:
        if line.startswith('#'):
            continue
        new = line.strip().split('\t')
        if new[2] != 'CDS':
            continue
        n = new[-1].split(';')
        for i, j in enumerate(n):
            if 'Parent=' in j:
                mindex = i
        g = n[mindex].split('.')
        t = g[0].replace('Parent=', '')
        if '_' in t:
            gene = t.split('_')[0]
        else:
            gene = t
        if gene not in mdict:
            mdict[gene] = []
        start1 = int(new[3])
        stop1 = int(new[4])
        mdict[gene].append((new[0], start1, stop1, new[6]))

    for i in sorted(mdict):
        k = ''
        for j in range(len(mdict[i])):
            k1 = f.sequence({'chr': mdict[i][j][0], 'start': mdict[i][j][1], 'stop': mdict[i][j][2], 'strand': mdict[i][j][3]})
            k += k1
        out1.write('>{0}-exon'.format(i) + '\n')
        out1.write(k + '\n')
    fh.close()
    out1.close()


def intron(fa, ann):
    f = Fasta(fa)
    fh = open(ann, 'r')
    out1 = open('intron.fa', 'w')
    mdict = {}
    ndict = {}
    for line in fh:
        if line.startswith('#'):
            continue
        new = line.strip().split('\t')
        if new[2] != 'CDS':
            continue
        n = new[-1].split(';')
        for i, j in enumerate(n):
            if 'Parent=' in j:
                mindex = i
        g = n[mindex].split('.')
        t = g[0].replace('Parent=', '')
        if '_' in t:
            gene = t.split('_')[0]
        else:
            gene = t
        if gene not in mdict:
            mdict[gene] = []
            ndict[gene] = [new[0], new[6]]
        start1 = int(new[3])
        stop1 = int(new[4])
        mdict[gene].append((start1, stop1))

    for i in sorted(mdict):
        k = ''
        total = len(mdict[i])
        for j in range(0, total - 1):
            start = mdict[i][j][1] + 1
            stop = mdict[i][j + 1][0] - 1
            k1 = f.sequence({'chr': ndict[i][0], 'start': start, 'stop': stop, 'strand': ndict[i][1]})
            k += k1
        out1.write('>{0}-intron'.format(i) + '\n')
        out1.write(k + '\n')
    fh.close()
    out1.close()
