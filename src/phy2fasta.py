#!/usr/bin/env python

__author__ = 'atamuri'

import sys
from Bio import SeqIO, AlignIO
from Bio.AlignIO import PhylipIO

if __name__ == "__main__":
    if len(sys.argv) == 1:
        alignment_in = AlignIO.read(sys.stdin, format='phylip-relaxed')
    else:
        alignment_in = AlignIO.read(sys.argv[1], format='phylip-relaxed')

    print(alignment_in.format('fasta'))

