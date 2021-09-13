#!/usr/bin/env python

__author__ = 'atamuri'

import sys
from Bio import SeqIO, AlignIO
from Bio.AlignIO import PhylipIO

def write_phylip_relaxed(handle, seqs):
    id_width = max((len(s.id.strip()) for s in seqs)) + 4
    alignment_out = AlignIO.MultipleSeqAlignment(seqs)
    writer = PhylipIO.SequentialPhylipWriter(handle)
    writer.write_alignment(alignment_out, id_width=id_width)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        alignment_in = SeqIO.parse(sys.stdin, format='fasta')
    else:
        alignment_in = SeqIO.parse(sys.argv[1], format='fasta')
    handle = sys.stdout
    write_phylip_relaxed(handle, list(alignment_in))
    handle.flush()


