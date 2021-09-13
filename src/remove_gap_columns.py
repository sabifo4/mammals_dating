#!/usr/bin/env python

__author__ = 'atamuri'

import sys
from Bio import SeqIO

def get_index_of_gaps(seq):
    index_of_gaps = set()
    for i in range(0, len(seq)):
        if seq[i] in ['-', 'N', '?']:
            index_of_gaps.add(i)
    return index_of_gaps

def get_all_gap_sites(seqs):
    index_of_all_gaps = set(range(0, len(seqs[0])))
    for seq in seqs:
        gap_sites = get_index_of_gaps(seq)
        if len(gap_sites) > 0:
            index_of_all_gaps = index_of_all_gaps.intersection(gap_sites)
    index_of_all_gaps = list(index_of_all_gaps)
    index_of_all_gaps.sort()
    return index_of_all_gaps

def remove_gaps(seqs, index_of_gaps):
    index_of_gaps.sort()
    for s in seqs:
        raw = s.seq
        for i in reversed(index_of_gaps):
            raw = raw[:i] + raw[i+1:]
        s.seq = raw

if __name__ == "__main__":
    if len(sys.argv) == 1:
        alignment = SeqIO.parse(sys.stdin, format='fasta')
    else:
        alignment = SeqIO.parse(sys.argv[1], format='fasta')

    alignment = list(alignment)
    all_gap_sites = get_all_gap_sites(alignment)
    if len(all_gap_sites) > 0:
        sys.stderr.write("Removing %s all gap column(s): %s\n" % (len(all_gap_sites), all_gap_sites))
        sys.stderr.flush()
        remove_gaps(alignment, all_gap_sites)

    SeqIO.write(alignment, sys.stdout, format='fasta')
    sys.stdout.flush()

