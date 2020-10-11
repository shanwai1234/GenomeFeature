#!/usr/bin/env python

import argparse
from GenomeFeature import GCTA
from GenomeFeature import surrounding
import argcomplete

__author__ = "Zhikai Liang"
__copyright__ = "Copyright 2019, GenomeFeature"
__license__ = "BSD 3.0"


def main():
    '''
    main code to calculate genome features for each gene
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--annotation', required=True,
                        help='This is the genome annotation file')
    parser.add_argument('-g', '--genome', required=True,
                        help='This is the genome assembly sequence')
    parser.add_argument('-d', '--fiveutr', required=True,
                        help='This is the sequence length for 5utr')
    parser.add_argument('-u', '--threeutr', required=True,
                        help='This is the sequence length for 3utr')
    parser.add_argument('-e', '--extension', required=True,
                        help='This is the extended sequence length for both upstream and downstream')
    parser.add_argument('-o', '--outfile', required=True,
                        help='This is the outfile with features for each gene')
    argcomplete.autocomplete(parser)
    args = parser.parse_args()
    utr5 = int(args.fiveutr)
    utr3 = int(args.threeutr)
    ext = int(args.extension)
    
    print 'Starting to generate fasta file for upstream regions of single genes ...\n'
    surrounding.upstream(args.genome, args.annotation, utr5, ext)
    print 'Starting to generate fasta file for downstream regions of single genes ...\n'
    surrounding.downstream(args.genome, args.annotation, utr3, ext)
    print 'Starting to generate fasta file for exon regions of single genes ...\n'
    surrounding.exon(args.genome, args.annotation)
    print 'Starting to generate fasta file for intron regions of single genes ...\n'
    surrounding.intron(args.genome, args.annotation)

    print 'Starting to calculate genome features in exon region ...\n'
    e1, e2 = GCTA.nucleotidefeature('exon.fa', 'Exon')
    tlist = ['Gene']
    tlist.extend(e1[1:])
    print 'Starting to calculate genome features in intron region ...\n'
    i1, i2 = GCTA.nucleotidefeature('intron.fa', 'Intron')
    tlist.extend(i1[1:])
    print 'Starting to calculate genome features in 5utr region ...\n'
    f1, f2 = GCTA.nucleotidefeature('5-UTR.fa', '5UTR')
    tlist.extend(f1[1:])
    print 'Starting to calculate genome features in 3utr region ...\n'
    t1, t2 = GCTA.nucleotidefeature('3-UTR.fa', '3UTR')
    tlist.extend(t1[1:])
    print 'Starting to calculate genome features in gene upstream region ...\n'
    u1, u2 = GCTA.nucleotidefeature('upstream.fa', 'Up')
    tlist.extend(u1[1:])
    print 'Starting to calculate genome features in gene downstream region ...\n'
    d1, d2 = GCTA.nucleotidefeature('downstream.fa', 'Down')
    tlist.extend(d1[1:])
    print 'Merging features ...\n'
    ndict = GCTA.merge(e2, i2, f2, t2, u2, d2)
    print 'Generating the final file with all features ...\n'
    mylen = len(tlist) - 1
    with open(args.outfile, 'w') as out:
        out.write(','.join(tlist) + '\n')
        for i in sorted(ndict):
            # only save genes with full sets of features
            if len(ndict[i]) != mylen:
                continue
            out.write(i + ',' + ','.join(ndict[i]) + '\n')


if __name__ == "__main__":
    main()
